/**
 * dsh-singlecell-workbench — Host 业务实现（唯一操作实现）
 *
 * 一条操作只写一次：浏览器 UI 走 HTTP 路由，Agent 走 singlecell 工具，
 * 两边的入口都调本文件里的 op* 函数（见 references/user-actions.md 的「一个操作两个调用方」）。
 *
 * 本文件按需用带 rev 的 URL 动态 import（见 index.js 的 loadImpl），
 * 因此不要在 lib 模块之间写静态相对 import；需要 engine.mjs 时用 loadEngine()。
 */

import { createReadStream, createWriteStream, statSync, unlinkSync } from 'node:fs'
import path from 'node:path'
import { randomUUID } from 'node:crypto'
import { pipeline } from 'node:stream/promises'
import { execFile } from 'node:child_process'

const PREFIX = '/singlecell'
const MAX_JSON_BODY = 1024 * 1024
const MAX_UPLOAD_BYTES = 8 * 1024 * 1024 * 1024 // 8 GiB，h5ad 可能很大
const ACTIVITY_LIMIT = 200
const STATUS_POLL = 50

/** 流程步骤定义：顺序即依赖，与 CONTRACT.md §2.2 一致。 */
export const STEPS = [
  { id: 'qc', title: '质控 QC', requires: [], params: { min_genes: 200, max_genes: null, max_pct_mt: 20 } },
  { id: 'filter', title: '过滤细胞', requires: ['qc'], params: { min_genes: 200, max_genes: null, min_counts: null, max_pct_mt: 20, min_cells: 3 } },
  { id: 'normalize', title: '标准化', requires: ['filter'], params: { target_sum: 10000 } },
  { id: 'hvg', title: '高变基因', requires: ['normalize'], params: { n_top_genes: 2000 } },
  { id: 'pca', title: 'PCA', requires: ['hvg'], params: { n_pcs: 30 } },
  { id: 'neighbors', title: '邻接图', requires: ['pca'], params: { n_neighbors: 15, metric: 'cosine' } },
  { id: 'cluster', title: '聚类', requires: ['neighbors'], params: { resolution: 1.0, method: 'louvain' } },
  { id: 'embed', title: '二维嵌入', requires: ['pca'], params: { layout: 'tsne', perplexity: 30 } },
  { id: 'markers', title: '标记基因', requires: ['cluster'], params: { top_n: 10, method: 't-test' } },
  { id: 'annotate', title: '细胞类型注释', requires: ['cluster'], params: { mapping: {} } },
]
export const STEP_IDS = STEPS.map((s) => s.id)
export const PLOT_KINDS = ['qc_violin', 'embedding', 'dotplot', 'heatmap', 'cluster_bar', 'pca_variance', 'qc_scatter']
const COLORS = ['cluster', 'obs', 'gene']

// ---------------------------------------------------------------------------
// HTTP 小工具
// ---------------------------------------------------------------------------

const LOOPBACK_HOSTS = new Set(['127.0.0.1', 'localhost', '[::1]', '::1', '0.0.0.0'])

export function hostIsLoopback(req) {
  const raw = String((req && req.headers && req.headers.host) || '')
  if (!raw) return true
  const name = raw.startsWith('[') ? raw.slice(0, raw.indexOf(']') + 1) : raw.split(':')[0]
  return LOOPBACK_HOSTS.has(name.toLowerCase())
}

function sendJson(res, status, value) {
  const body = Buffer.from(JSON.stringify(value), 'utf8')
  res.writeHead(status, {
    'content-type': 'application/json; charset=utf-8',
    'cache-control': 'no-store',
    'content-length': String(body.length),
  })
  res.end(body)
}

function fail(code, message, status = 400, details = null) {
  const body = { ok: false, code, error: message }
  if (details && typeof details === 'object') Object.assign(body, details)
  return { status, body }
}

function toFailure(error) {
  const code = (error && error.code) || 'internal'
  const message = String((error && error.message) || error)
  const details = error && error.details && typeof error.details === 'object' ? error.details : null
  if (code === 'ambiguous') return fail(code, message, 400, details)
  if (code === 'bad_request' || code === 'not_found' || code === 'prerequisite' || code === 'unsupported') {
    return fail(code, message, 400, details)
  }
  if (code === 'timeout') return fail(code, message, 504, details)
  if (code === 'no_python' || code === 'engine_exited' || code === 'internal') return fail(code, message, 500, details)
  if (code === 'busy') return fail(code, message, 429, details)
  return fail(code, message, 500, details)
}

function readBody(req, limit) {
  return new Promise((resolve, reject) => {
    const chunks = []
    let size = 0
    req.on('data', (chunk) => {
      size += chunk.length
      if (size > limit) {
        reject(Object.assign(new Error(`请求体超过上限（${Math.round(limit / 1024)} KB）`), { code: 'bad_request' }))
        req.destroy()
        return
      }
      chunks.push(chunk)
    })
    req.on('end', () => resolve(Buffer.concat(chunks)))
    req.on('error', (error) => reject(Object.assign(error, { code: 'bad_request' })))
  })
}

async function readJson(req) {
  const buf = await readBody(req, MAX_JSON_BODY)
  if (!buf.length) return {}
  try {
    const value = JSON.parse(buf.toString('utf8'))
    if (value === null || typeof value !== 'object' || Array.isArray(value)) {
      throw new Error('请求体必须是 JSON 对象')
    }
    return value
  } catch (error) {
    throw Object.assign(new Error(`请求体不是合法 JSON：${error && error.message ? error.message : error}`), { code: 'bad_request' })
  }
}

function sanitizeName(name, fallback = 'data.h5ad') {
  const base = String(name || '').split(/[\\/]/).pop() || fallback
  const cleaned = base.replace(/[^\w.\-]+/g, '_').slice(0, 120)
  return cleaned || fallback
}

// ---------------------------------------------------------------------------
// runtime / 状态
// ---------------------------------------------------------------------------

export function stateOf(runtime) {
  if (!runtime.sc) {
    runtime.sc = {
      active: null,
      datasets: new Map(),
      activity: [],
      figures: new Map(),
      suggestCache: { key: '', value: [], at: 0 },
      lastRefresh: 0,
      bootAt: Date.now(),
    }
  }
  return runtime.sc
}

export function logActivity(runtime, entry) {
  const sc = stateOf(runtime)
  sc.activity.push({
    at: new Date().toISOString(),
    actor: entry.actor === 'agent' ? 'agent' : entry.actor === 'system' ? 'system' : 'user',
    action: String(entry.action || 'unknown'),
    dataset: entry.dataset || null,
    detail: String(entry.detail || '').slice(0, 400),
    ok: entry.ok !== false,
    ms: typeof entry.ms === 'number' ? Math.round(entry.ms) : null,
  })
  if (sc.activity.length > ACTIVITY_LIMIT) sc.activity.splice(0, sc.activity.length - ACTIVITY_LIMIT)
}

export function activityOf(runtime, limit = STATUS_POLL) {
  const sc = stateOf(runtime)
  return sc.activity.slice(-limit).reverse()
}

async function loadEngineMod() {
  return import(new URL(`./engine.mjs?rev=${Date.now()}-${Math.random().toString(36).slice(2)}`, import.meta.url).href)
}

/** 拿到引擎实例；必要时探测 Python 并启动。失败结果缓存 60s，避免 UI 轮询时反复探测。 */
async function ensureEngine(runtime, { allowStart = true } = {}) {
  if (runtime.engine) {
    if (runtime.engine.running) return runtime.engine
    if (!allowStart) return runtime.engine
  } else {
    const mod = await loadEngineMod()
    const config = runtime.config && typeof runtime.config === 'object' ? runtime.config : {}
    runtime.engine = mod.createEngine({
      python: typeof config.python === 'string' && config.python.trim() ? config.python.trim() : undefined,
      workdir: typeof config.workdir === 'string' && config.workdir.trim() ? config.workdir.trim() : undefined,
      timeout: Number(config.timeoutMs) > 0 ? Number(config.timeoutMs) : undefined,
      config,
    })
    if (typeof runtime.ctx?.effect === 'function') {
      const engine = runtime.engine
      runtime.ctx.effect(() => () => { engine.dispose().catch(() => {}) }, 'dsh-singlecell-workbench: engine')
    }
  }
  if (!allowStart) return runtime.engine
  const now = Date.now()
  if (runtime.engineFail && now - runtime.engineFail.at < 60_000) throw runtime.engineFail.error
  try {
    await runtime.engine.start()
    runtime.engineFail = null
  } catch (error) {
    runtime.engineFail = { at: now, error }
    throw error
  }
  return runtime.engine
}

/** 引擎概要：永远不抛异常，供 /health 与 /state 展示。 */
export async function engineStatus(runtime) {
  let engine = runtime.engine
  if (!engine) {
    try {
      await ensureEngine(runtime, { allowStart: false })
      engine = runtime.engine
    } catch (error) {
      return { ready: false, mode: 'none', python: null, error: String(error && error.message ? error.message : error) }
    }
  }
  const status = engine.status()
  const mode = status.mode === 'unknown' ? (status.error ? 'none' : 'builtin') : status.mode
  return {
    ready: Boolean(status.ready),
    mode: status.ready ? mode : status.error ? 'none' : 'unknown',
    python: status.python,
    version: status.version,
    packages: status.packages,
    layouts: status.layouts,
    warnings: status.warnings || [],
    error: status.error || null,
    progress: status.lastProgress ? status.lastProgress.message : null,
  }
}

function syncDataset(runtime, dataset, { activate = false } = {}) {
  if (!dataset || !dataset.id) return null
  const sc = stateOf(runtime)
  const previous = sc.datasets.get(dataset.id) || {}
  const merged = { ...previous, ...dataset }
  sc.datasets.set(dataset.id, merged)
  if (activate || !sc.active) sc.active = dataset.id
  return merged
}

function datasetList(runtime) {
  const sc = stateOf(runtime)
  return Array.from(sc.datasets.values())
}

async function refreshDatasets(runtime, { force = false } = {}) {
  const sc = stateOf(runtime)
  const now = Date.now()
  if (!force && now - sc.lastRefresh < 1200) return
  sc.lastRefresh = now
  let engine
  try {
    engine = await ensureEngine(runtime)
  } catch {
    // 引擎不可用时清掉内存中已失效的数据集（worker 里的数据已经没了）
    sc.datasets.clear()
    sc.active = null
    return
  }
  try {
    const result = await engine.request('list', {}, { timeout: 15_000 })
    const datasets = Array.isArray(result?.datasets) ? result.datasets : []
    const ids = new Set(datasets.map((d) => d.id))
    for (const d of datasets) syncDataset(runtime, d)
    for (const id of Array.from(sc.datasets.keys())) if (!ids.has(id)) sc.datasets.delete(id)
    if (sc.active && !ids.has(sc.active)) sc.active = datasets.length ? datasets[0].id : null
  } catch {
    // 保留上次缓存，不打扰用户
  }
}

function resolveDatasetId(runtime, dataset) {
  const sc = stateOf(runtime)
  if (dataset && sc.datasets.has(dataset)) return dataset
  if (sc.active) return sc.active
  const first = datasetList(runtime)[0]
  if (first) return first.id
  throw Object.assign(new Error('还没有载入任何数据集，请先载入示例数据或打开 .h5ad'), { code: 'bad_request' })
}

// ---------------------------------------------------------------------------
// 建议下一步（确定性，不调用模型）
// ---------------------------------------------------------------------------

async function suggestionsFor(runtime, datasetId) {
  const sc = stateOf(runtime)
  const ds = datasetId ? sc.datasets.get(datasetId) : null
  if (!ds) {
    return [
      { id: 'demo', title: '载入示例数据', detail: '内置 2400 细胞 × 1200 基因的合成数据，可直接跑完整流程', action: 'demo' },
      { id: 'open', title: '打开自己的 h5ad', detail: '在顶栏输入 .h5ad 绝对路径，或直接上传文件', action: 'open' },
    ]
  }
  const key = `${ds.id}:${(ds.steps || []).join(',')}:${ds.n_obs}x${ds.n_vars}`
  if (sc.suggestCache.key === key && Date.now() - sc.suggestCache.at < 5000) return sc.suggestCache.value
  let value = []
  try {
    const engine = await ensureEngine(runtime)
    const result = await engine.request('suggest', { dataset: ds.id }, { timeout: 20_000 })
    value = Array.isArray(result?.suggestions) ? result.suggestions : []
  } catch {
    value = []
  }
  sc.suggestCache = { key, value, at: Date.now() }
  return value
}

/** 预设流程（skills 的本地版）：引擎提供定义，Host 负责顺序执行。 */
async function recipesFor(runtime) {
  const sc = stateOf(runtime)
  if (sc.recipeCache && Date.now() - sc.recipeCache.at < 60_000) return sc.recipeCache.value
  let value = []
  try {
    const engine = await ensureEngine(runtime)
    const result = await engine.request('recipe', {}, { timeout: 20_000 })
    value = Array.isArray(result?.recipes) ? result.recipes : []
  } catch {
    value = []
  }
  sc.recipeCache = { value, at: Date.now() }
  return value
}

export function recipeProgressOf(runtime) {
  const sc = stateOf(runtime)
  return sc.recipeProgress || null
}

// ---------------------------------------------------------------------------
// 操作（UI 与 Agent 工具共用）
// ---------------------------------------------------------------------------

function actorOf(value, fallback = 'user') {
  return value === 'agent' || value === 'system' ? value : fallback
}

export async function opHealth(runtime) {
  const engine = await engineStatus(runtime)
  if (!engine.ready && !engine.error) {
    try { await ensureEngine(runtime) } catch (error) {
      engine.error = String(error && error.message ? error.message : error)
      engine.mode = 'none'
    }
    const next = await engineStatus(runtime)
    Object.assign(engine, next)
  }
  return { ok: true, engine }
}

export async function opState(runtime) {
  await refreshDatasets(runtime)
  const sc = stateOf(runtime)
  const engine = await engineStatus(runtime)
  const active = sc.active
  const suggestions = await suggestionsFor(runtime, active)
  const recipes = await recipesFor(runtime)
  return {
    ok: true,
    engine,
    active,
    datasets: datasetList(runtime).map((d) => ({
      id: d.id,
      name: d.name || d.id,
      path: d.path || null,
      n_obs: d.n_obs,
      n_vars: d.n_vars,
      steps: Array.isArray(d.steps) ? d.steps : [],
      current: d.current || null,
      warnings: Array.isArray(d.warnings) ? d.warnings : [],
      obs_columns: Array.isArray(d.obs_columns) ? d.obs_columns : [],
      samples: Array.isArray(d.samples) ? d.samples : [],
      groups: Array.isArray(d.groups) ? d.groups : [],
      group_column: d.group_column || null,
      history: Array.isArray(d.history) ? d.history.slice(-20) : [],
      source_format: d.source_format || d.source || null,
      notes: Array.isArray(d.notes) ? d.notes : [],
      analyses_run: Array.isArray(d.analyses_run) ? d.analyses_run.slice(-20) : [],
    })),
    suggestions,
    recipes,
    recipeProgress: recipeProgressOf(runtime),
    activity: activityOf(runtime, 60),
  }
}

export async function opDemo(runtime, { actor } = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  const result = await engine.request('open', { demo: true }, { timeout: 120_000 })
  const dataset = syncDataset(runtime, result.dataset, { activate: true })
  logActivity(runtime, {
    actor: actorOf(actor), action: 'open:demo', dataset: dataset.id, ms: Date.now() - started,
    detail: `${dataset.name || 'demo'} · ${dataset.n_obs} 细胞 × ${dataset.n_vars} 基因`,
  })
  return { ok: true, dataset, datasets: datasetList(runtime) }
}

export async function opOpen(runtime, { path: filePath, url, name, format, member, object, transpose, paths, mergeDir, samples, actor } = {}) {
  const started = Date.now()
  const target = String(filePath || '').trim()
  const remote = String(url || '').trim()
  const hasPaths = Array.isArray(paths) && paths.length > 0
  const hasSamples = Array.isArray(samples) && samples.length > 0
  if (!target && !remote && !hasPaths && !hasSamples) {
    throw Object.assign(new Error('请提供本地路径、URL、paths 多文件合并，或 samples 样本表（扫描根目录后勾选样本得来）'), { code: 'bad_request' })
  }
  const engine = await ensureEngine(runtime)
  const args = {}
  if (target) args.path = target
  if (remote) args.url = remote
  if (name) args.name = name
  if (format) args.format = format
  if (member) args.member = member
  if (object) args.object = object
  if (transpose !== undefined && transpose !== null) args.transpose = Boolean(transpose)
  if (Array.isArray(paths) && paths.length) args.paths = paths.map(String)
  if (mergeDir) args.merge_dir = true
  if (hasSamples) {
    args.samples = samples.map((item) => ({
      path: String(item.path || item.source || ''),
      name: item.name ? String(item.name) : undefined,
      group: item.group ? String(item.group) : undefined,
      trio_prefix: item.trio_prefix == null ? undefined : String(item.trio_prefix),
      format: item.format ? String(item.format) : undefined,
    })).filter((item) => item.path)
  }
  let result
  try {
    result = await engine.request('open', args, { timeout: 1_800_000 })
  } catch (error) {
    logActivity(runtime, {
      actor: actorOf(actor), action: 'open', detail: `${remote || target} → ${error.message}`, ok: false, ms: Date.now() - started,
    })
    throw error
  }
  const dataset = syncDataset(runtime, result.dataset, { activate: true })
  logActivity(runtime, {
    actor: actorOf(actor), action: 'open', dataset: dataset.id, ms: Date.now() - started,
    detail: `${dataset.name} · ${dataset.source_format || 'auto'} · ${dataset.n_obs} 细胞 × ${dataset.n_vars} 基因`,
  })
  return { ok: true, dataset, datasets: datasetList(runtime) }
}

const PICKER_TIMEOUT_MS = 15 * 60 * 1000

/** 构造各平台的原生选择框命令（不要拼接用户可控的 shell，参数以数组传递）。 */
function pickerCommand(kind, prompt) {
  const wantDir = String(kind || 'dir') !== 'files'
  const text = String(prompt || (wantDir ? '选择数据根目录' : '选择数据文件')).replace(/[\r\n"]/g, ' ').slice(0, 200)
  if (process.platform === 'darwin') {
    const script = wantDir
      ? `POSIX path of (choose folder with prompt "${text}")`
      : `set out to {}\nrepeat with f in (choose file with prompt "${text}" with multiple selections allowed)\nset end of out to POSIX path of f\nend repeat\nreturn out`
    return { command: '/usr/bin/osascript', args: ['-e', script], kind: wantDir ? 'dir' : 'files' }
  }
  if (process.platform === 'win32') {
    const script = wantDir
      ? `Add-Type -AssemblyName System.Windows.Forms; $d=New-Object System.Windows.Forms.FolderBrowserDialog; $d.Description='${text}'; if($d.ShowDialog() -eq 'OK'){Write-Output $d.SelectedPath}`
      : `Add-Type -AssemblyName System.Windows.Forms; $d=New-Object System.Windows.Forms.OpenFileDialog; $d.Multiselect=$true; $d.Title='${text}'; if($d.ShowDialog() -eq 'OK'){$d.FileNames}`
    return { command: 'powershell.exe', args: ['-NoProfile', '-STA', '-Command', script], kind: wantDir ? 'dir' : 'files' }
  }
  const args = wantDir
    ? ['--file-selection', '--directory', '--title=' + text]
    : ['--file-selection', '--multiple', '--separator=\n', '--title=' + text]
  return { command: 'zenity', args, kind: wantDir ? 'dir' : 'files' }
}

/** 打开本机原生选择框，返回绝对路径（浏览器无法提供绝对路径，必须由 Host 弹框）。 */
export async function opPickPath(runtime, { kind, prompt, dryRun } = {}) {
  const spec = pickerCommand(kind, prompt)
  if (dryRun) return { ok: true, dry_run: true, platform: process.platform, command: spec.command, args: spec.args, kind: spec.kind }
  const started = Date.now()
  const result = await new Promise((resolve) => {
    execFile(spec.command, spec.args, { timeout: PICKER_TIMEOUT_MS, windowsHide: true }, (error, stdout) => {
      if (error) {
        const code = typeof error.code === 'number' ? error.code : null
        if (code === 1 || error.killed || error.signal) {
          return resolve({ ok: true, cancelled: true, paths: [] })
        }
        return resolve({ ok: false, cancelled: false, paths: [], error: `无法打开系统选择框：${error.message}` })
      }
      const paths = String(stdout || '')
        .split(/\r?\n/)
        .map((line) => line.trim())
        .filter(Boolean)
      resolve({ ok: true, cancelled: paths.length === 0, paths })
    })
  })
  logActivity(runtime, {
    actor: actorOf('user'), action: 'pick',
    detail: result.ok ? (result.cancelled ? '已取消选择' : `${spec.kind === 'dir' ? '目录' : '文件'}：${result.paths[0] || ''}`) : result.error,
    ok: Boolean(result.ok), ms: Date.now() - started,
  })
  return result
}

/** 扫描 GEO 数据根目录，自动识别样本（10X 目录 / GEO 三联体 / h5 / h5ad / rds / RData / 矩阵）。 */
/** 提取小亚群：按簇子集化生成新数据集（Step4 二次降维用）。 */
export async function opSubset(runtime, { dataset, clusters, column, name, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const list = (Array.isArray(clusters) ? clusters : [clusters])
    .map((item) => String(item == null ? '' : item)).filter(Boolean)
  if (!list.length) throw Object.assign(new Error('请提供要提取的簇（clusters: ["3","7"]）'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const body = { dataset: id, clusters: list, column: column || undefined, name: name || undefined }
  const result = await engine.request('subset', body, { timeout: 1_800_000 })
  const datasetResult = result && result.dataset ? syncDataset(runtime, result.dataset, { activate: true }) : null
  logActivity(runtime, {
    actor: actorOf(actor), action: 'subset', dataset: id,
    detail: `提取簇 ${list.join('、')} → ${result && result.subset_cells || 0} 个细胞`,
  })
  return { ok: true, dataset: datasetResult, parent: result && result.parent, subset_cells: result && result.subset_cells, kept_clusters: list }
}

export async function opScan(runtime, { path, recursive, actor } = {}) {
  const root = String(path || '').trim()
  if (!root) throw Object.assign(new Error('请提供数据根目录路径（GEO 解压后的文件夹）'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const result = await engine.request('scan', { path: root, recursive: Boolean(recursive) }, { timeout: 300_000 })
  logActivity(runtime, {
    actor: actorOf(actor), action: 'scan', detail: `${root} → ${(result.samples || []).length} 个样本`, ok: true,
  })
  return { ok: true, ...result }
}

/** GEO 编号查询：标题 + 每个 GSM 的说明，辅助样本分组。 */
export async function opGeoLookup(runtime, { gse } = {}) {
  const value = String(gse || '').trim().toUpperCase()
  if (!value) throw Object.assign(new Error('请提供 GSE 编号（如 GSE246662）'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const result = await engine.request('geo_lookup', { gse: value }, { timeout: 60_000 })
  return { ok: true, ...result }
}

/** 只探测格式/维度，不加载矩阵（大文件/压缩包歧义时用）。 */
export async function opInspect(runtime, { path: filePath, url, format, member, object } = {}) {
  const target = String(filePath || '').trim()
  const remote = String(url || '').trim()
  if (!target && !remote) throw Object.assign(new Error('请提供本地路径或 URL'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const args = {}
  if (target) args.path = target
  if (remote) args.url = remote
  if (format) args.format = format
  if (member) args.member = member
  if (object) args.object = object
  const result = await engine.request('inspect', args, { timeout: 300_000 })
  return { ok: true, ...result }
}

/** 支持的导入格式清单（含 R 桥是否可用）。 */
export async function opFormats(runtime) {
  const engine = await ensureEngine(runtime)
  const result = await engine.request('formats', {}, { timeout: 30_000 })
  return { ok: true, formats: Array.isArray(result?.formats) ? result.formats : [] }
}

/** 高级分析目录（native + external，含可用性与参数描述）。 */
export async function opAnalyses(runtime, { dataset } = {}) {
  const engine = await ensureEngine(runtime)
  const sc = stateOf(runtime)
  const id = dataset || sc.active || undefined
  const result = await engine.request('analyses', id ? { dataset: id } : {}, { timeout: 120_000 })
  const analyses = Array.isArray(result?.analyses) ? result.analyses : []
  return {
    ok: true,
    dataset: result?.dataset || id || null,
    groups: Array.isArray(result?.groups) ? result.groups : [],
    columns: Array.isArray(result?.columns) ? result.columns : [],
    group_keys: Array.isArray(result?.group_keys) ? result.group_keys : [],
    diagnostics: result?.diagnostics || {},
    analyses: analyses.map((item) => ({
      ...item,
      available: Boolean(item.native && item.available),
      reason: item.reason || null,
    })),
  }
}

/** 跑一个高级分析；图注册进 figure 表，obs 列由引擎写回数据集。 */
export async function opAnalysis(runtime, { dataset, id: analysisId, params, actor } = {}) {
  const target = resolveDatasetId(runtime, dataset)
  const analysis = String(analysisId || '').trim()
  if (!analysis) throw Object.assign(new Error('请提供分析 id（先看 GET /singlecell/analyses）'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const safeParams = params && typeof params === 'object' ? params : {}
  let result
  try {
    // 官方工具（ext_*）可能是长任务：单步最多 2 小时
    const analysisTimeout = analysis.startsWith('ext_') ? 7_200_000 : 1_800_000
    result = await engine.request('analysis', { dataset: target, id: analysis, params: safeParams }, { timeout: analysisTimeout })
  } catch (error) {
    logActivity(runtime, {
      actor: actorOf(actor), action: `analysis:${analysis}`, dataset: target, ok: false,
      detail: error.message, ms: Date.now() - started,
    })
    throw error
  }
  const sc = stateOf(runtime)
  const figures = (Array.isArray(result.figures) ? result.figures : []).map((figure) => {
    const file = String(figure.file || '')
    const figureId = path.basename(file)
    if (!figureId || !file) return figure
    sc.figures.set(figureId, { file, width: figure.width, height: figure.height, kind: 'analysis', dataset: target, at: Date.now() })
    return { ...figure, id: figureId, url: `${PREFIX}/figure?id=${encodeURIComponent(figureId)}` }
  })
  const bullets = Array.isArray(result.summary?.bullets) ? result.summary.bullets : []
  logActivity(runtime, {
    actor: actorOf(actor), action: `analysis:${analysis}`, dataset: target, ms: Date.now() - started,
    detail: `${result.title || analysis}${bullets.length ? '：' + bullets.slice(0, 2).join('；') : ''}`,
  })
  return { ok: true, ...result, figures, dataset: target }
}

export async function opUpload(runtime, { filename, buffer, actor } = {}) {
  const started = Date.now()
  if (!buffer || !buffer.length) throw Object.assign(new Error('上传内容为空'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const safe = sanitizeName(filename)
  const dest = path.join(engine.dirs.uploads, `${Date.now().toString(36)}-${randomUUID().slice(0, 8)}-${safe}`)
  await new Promise((resolve, reject) => {
    const out = createWriteStream(dest)
    out.on('error', reject)
    out.on('finish', resolve)
    out.end(buffer)
  })
  try {
    const result = await engine.request('open', { path: dest, name: safe }, { timeout: 180_000 })
    const dataset = syncDataset(runtime, result.dataset, { activate: true })
    logActivity(runtime, {
      actor: actorOf(actor), action: 'upload', dataset: dataset.id, ms: Date.now() - started,
      detail: `${safe} · ${dataset.n_obs} 细胞 × ${dataset.n_vars} 基因`,
    })
    return { ok: true, dataset, datasets: datasetList(runtime) }
  } catch (error) {
    logActivity(runtime, { actor: actorOf(actor), action: 'upload', detail: `${safe} → ${error.message}`, ok: false, ms: Date.now() - started })
    throw error
  }
}

export async function opSelect(runtime, { dataset } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  stateOf(runtime).active = id
  return { ok: true, active: id, datasets: datasetList(runtime) }
}

export async function opClose(runtime, { dataset, actor } = {}) {
  const sc = stateOf(runtime)
  const id = String(dataset || sc.active || '')
  if (!id) throw Object.assign(new Error('没有可关闭的数据集'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  try { await engine.request('close', { dataset: id }, { timeout: 30_000 }) } catch { /* 引擎侧已无该数据集也无妨 */ }
  sc.datasets.delete(id)
  if (sc.active === id) sc.active = datasetList(runtime)[0]?.id || null
  logActivity(runtime, { actor: actorOf(actor), action: 'close', dataset: id, detail: '已关闭数据集' })
  return { ok: true, active: sc.active, datasets: datasetList(runtime) }
}

export async function opPreview(runtime, { dataset } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const engine = await ensureEngine(runtime)
  const result = await engine.request('preview', { dataset: id }, { timeout: 30_000 })
  return { ok: true, dataset: id, ...result }
}

export async function opRun(runtime, { dataset, step, params, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const stepId = String(step || '')
  if (!STEP_IDS.includes(stepId)) {
    throw Object.assign(new Error(`未知步骤 "${step}"，可用：${STEP_IDS.join('/')}`), { code: 'bad_request' })
  }
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const timeout = stepId === 'embed' || stepId === 'cluster' || stepId === 'markers' ? 300_000 : 180_000
  let result
  try {
    result = await engine.request('run', { dataset: id, step: stepId, params: params || {} }, { timeout })
  } catch (error) {
    logActivity(runtime, {
      actor: actorOf(actor), action: `run:${stepId}`, dataset: id, ok: false,
      detail: error.message, ms: Date.now() - started,
    })
    throw error
  }
  const sc = stateOf(runtime)
  const ds = sc.datasets.get(id)
  if (ds) {
    ds.steps = Array.isArray(result.steps) ? result.steps : ds.steps
    ds.current = result.current ?? stepId
    if (Array.isArray(result.history)) ds.history = result.history
    if (result.dataset) {
      ds.n_obs = result.dataset.n_obs ?? ds.n_obs
      ds.n_vars = result.dataset.n_vars ?? ds.n_vars
    }
    sc.datasets.set(id, ds)
  }
  const bullets = Array.isArray(result.summary?.bullets) ? result.summary.bullets : []
  logActivity(runtime, {
    actor: actorOf(actor), action: `run:${stepId}`, dataset: id, ms: Date.now() - started,
    detail: bullets.join('；') || result.summary?.title || stepId,
  })
  const suggestions = await suggestionsFor(runtime, id)
  return { ok: true, dataset: id, step: stepId, summary: result.summary || {}, steps: result.steps || [], current: result.current || stepId, suggestions }
}

export async function opEmbedding(runtime, { dataset, color_by, max_points } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const engine = await ensureEngine(runtime)
  const colorBy = normalizeColorBy(color_by)
  const result = await engine.request(
    'embedding',
    { dataset: id, color_by: colorBy, max_points: Number(max_points) > 0 ? Number(max_points) : 12000 },
    { timeout: 120_000 },
  )
  return { ok: true, dataset: id, ...result }
}

function normalizeColorBy(colorBy) {
  if (!colorBy || typeof colorBy !== 'object') return { kind: 'cluster', key: 'clusters' }
  const kind = COLORS.includes(colorBy.kind) ? colorBy.kind : 'cluster'
  const key = String(colorBy.key || (kind === 'cluster' ? 'clusters' : kind === 'gene' ? '' : 'n_genes'))
  if (kind === 'gene' && !key) throw Object.assign(new Error('按基因着色需要提供基因名'), { code: 'bad_request' })
  return { kind, key }
}

export async function opPlot(runtime, { dataset, kind, params, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const plotKind = String(kind || '')
  if (!PLOT_KINDS.includes(plotKind)) {
    throw Object.assign(new Error(`未知图类型 "${kind}"，可用：${PLOT_KINDS.join('/')}`), { code: 'bad_request' })
  }
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const safeParams = params && typeof params === 'object' ? params : {}
  const result = await engine.request('plot', { dataset: id, kind: plotKind, params: safeParams }, { timeout: 180_000 })
  const file = String(result.file || '')
  const figureId = path.basename(file)
  const sc = stateOf(runtime)
  sc.figures.set(figureId, { file, width: result.width, height: result.height, kind: plotKind, dataset: id, at: Date.now() })
  logActivity(runtime, {
    actor: actorOf(actor), action: `plot:${plotKind}`, dataset: id, ms: Date.now() - started,
    detail: `${result.width}×${result.height}px`,
  })
  return {
    ok: true,
    dataset: id,
    figure: {
      id: figureId,
      url: `${PREFIX}/figure?id=${encodeURIComponent(figureId)}`,
      width: result.width,
      height: result.height,
      kind: plotKind,
      file,
    },
    params: safeParams,
  }
}

export async function opMarkers(runtime, { dataset, top_n, method, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const result = await engine.request(
    'markers',
    { dataset: id, groupby: 'clusters', top_n: Number(top_n) > 0 ? Math.min(Number(top_n), 50) : 10, method: method || 't-test' },
    { timeout: 180_000 },
  )
  const groups = Array.isArray(result.groups) ? result.groups : []
  logActivity(runtime, {
    actor: actorOf(actor), action: 'markers', dataset: id, ms: Date.now() - started,
    detail: `${groups.length} 个簇，每簇 ${groups[0]?.genes?.length || 0} 个基因`,
  })
  return { ok: true, dataset: id, groups }
}

/** 确定性自检（reviewer 的轻量版）：不调用模型，检查流程与数据是否可站得住。 */
async function callMeta(runtime, cmd, args) {
  const engine = await ensureEngine(runtime)
  return engine.request(cmd, args, { timeout: 120_000 })
}

export async function opObs(runtime, { dataset } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const result = await callMeta(runtime, 'obs_list', { dataset: id })
  return { ok: true, dataset: id, ...result }
}

export async function opSamplesUpdate(runtime, { dataset, samples, groups, sampleColumn, groupColumn, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const result = await callMeta(runtime, 'samples_update', {
    dataset: id, samples: samples || [], groups: groups,
    sample_column: sampleColumn || 'sample', group_column: groupColumn || 'group',
  })
  logActivity(runtime, { actor: actorOf(actor), action: 'samples:update', dataset: id,
    detail: `${(samples || []).length} 个样本；分组列=${result.group_column || 'group'}` })
  return { ok: true, dataset: id, ...result }
}

/** 按样本移除细胞（真正剔除，重置下游结果）；samples 里放样本名/id。 */
export async function opSamplesRemove(runtime, { dataset, samples, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const list = (Array.isArray(samples) ? samples : [samples])
    .map((item) => String(item == null ? '' : (typeof item === 'object' ? (item.name || item.id || '') : item)))
    .filter(Boolean)
  if (!list.length) throw Object.assign(new Error('请提供要移除的样本名（samples）'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const result = await engine.request('samples_remove', { dataset: id, samples: list }, { timeout: 1_200_000 })
  const datasetResult = result && result.dataset ? syncDataset(runtime, result.dataset, { activate: true }) : null
  logActivity(runtime, {
    actor: actorOf(actor), action: 'samples:remove', dataset: id,
    detail: `移除样本 ${list.join('、')}（${(result && result.removed_cells) || 0} 个细胞）`,
  })
  return {
    ok: true,
    dataset: datasetResult,
    removed_cells: (result && result.removed_cells) || 0,
    removed_samples: (result && result.removed_samples) || list,
  }
}

export async function opSamplesAutoGroup(runtime, { dataset, pattern, sampleColumn } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const result = await callMeta(runtime, 'samples_auto_group', { dataset: id, pattern, sample_column: sampleColumn || 'sample' })
  return { ok: true, dataset: id, ...result }
}

export async function opSamplesImportMap(runtime, { dataset, path: mapPath, sampleColumn, groupColumn } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const result = await callMeta(runtime, 'samples_import_map', {
    dataset: id, path: mapPath, sample_column: sampleColumn, group_column: groupColumn,
  })
  return { ok: true, dataset: id, ...result }
}

export async function opObsJoin(runtime, { dataset, path: metaPath, key } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const result = await callMeta(runtime, 'obs_join', { dataset: id, path: metaPath, key })
  return { ok: true, dataset: id, ...result }
}

export async function opObsDerive(runtime, { dataset, name, source, pattern, group, value } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const result = await callMeta(runtime, 'obs_derive', { dataset: id, name, source, pattern, group, value })
  return { ok: true, dataset: id, ...result }
}

export async function opObsSetGroup(runtime, { dataset, column } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const result = await callMeta(runtime, 'obs_set_group', { dataset: id, column })
  return { ok: true, dataset: id, ...result }
}

export async function opAudit(runtime, { dataset, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const result = await engine.request('audit', { dataset: id }, { timeout: 60_000 })
  const checks = Array.isArray(result.checks) ? result.checks : []
  const problems = checks.filter((c) => c.level !== 'ok')
  logActivity(runtime, {
    actor: actorOf(actor), action: 'audit', dataset: id, ms: Date.now() - started,
    detail: `score=${result.score ?? '?'}，${problems.length} 项提示`,
  })
  return { ok: true, dataset: id, score: result.score ?? 0, level: result.level || 'ok', checks }
}

/** 参考 marker 面板的自动注释建议（只给建议，不写 obs）。 */
export async function opAutoannotate(runtime, { dataset, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const result = await engine.request('autoannotate', { dataset: id, groupby: 'clusters' }, { timeout: 120_000 })
  const suggestions = Array.isArray(result.suggestions) ? result.suggestions : []
  logActivity(runtime, {
    actor: actorOf(actor), action: 'autoannotate', dataset: id, ms: Date.now() - started,
    detail: suggestions.slice(0, 3).map((s) => `${s.cluster}→${s.label}`).join('，') || '无建议',
  })
  return { ok: true, dataset: id, panel: result.panel || 'builtin-marker-panel-v1', suggestions }
}

function normalizeResolutions(input) {
  const raw = Array.isArray(input) ? input : String(input == null ? '' : input).split(/[,，\s]+/)
  const nums = raw.map((v) => Number(v)).filter((v) => Number.isFinite(v) && v > 0 && v <= 5)
  if (!nums.length) {
    throw Object.assign(new Error('请提供 0.1–5 之间的分辨率，例如 "0.4,0.8,1.2"'), { code: 'bad_request' })
  }
  return Array.from(new Set(nums)).slice(0, 8)
}

/** 参数对比：在既有 kNN 图上扫描不同 resolution，看簇数与 ARI（不改当前聚类）。 */
export async function opResolutionScan(runtime, { dataset, resolutions, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const list = normalizeResolutions(resolutions)
  const result = await engine.request('resolution_scan', { dataset: id, resolutions: list }, { timeout: 240_000 })
  logActivity(runtime, {
    actor: actorOf(actor), action: 'resolution_scan', dataset: id, ms: Date.now() - started,
    detail: list.join('/'),
  })
  return { ok: true, dataset: id, current: result.current || null, results: Array.isArray(result.results) ? result.results : [] }
}

/** 一键跑预设流程（skills 的本地版）：逐步调 opRun，失败即停并保留已完成步骤。 */
export async function opRecipe(runtime, { dataset, id: recipeId, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const wanted = String(recipeId || 'standard')
  const recipes = await recipesFor(runtime)
  const recipe = recipes.find((r) => r.id === wanted)
  if (!recipe) {
    throw Object.assign(
      new Error(`未知预设流程 "${wanted}"，可用：${recipes.map((r) => r.id).join('/') || '（引擎未返回预设）'}`),
      { code: 'bad_request' },
    )
  }
  const steps = Array.isArray(recipe.steps) ? recipe.steps : []
  if (!steps.length) throw Object.assign(new Error(`预设 "${wanted}" 没有步骤`), { code: 'bad_request' })
  const sc = stateOf(runtime)
  const done = []
  try {
    for (let i = 0; i < steps.length; i += 1) {
      const current = steps[i]
      sc.recipeProgress = {
        recipe: wanted, title: recipe.title || wanted,
        index: i + 1, total: steps.length, step: current.step, at: new Date().toISOString(),
      }
      const result = await opRun(runtime, { dataset: id, step: current.step, params: current.params || {}, actor })
      done.push({ step: current.step, ok: true, title: result.summary?.title || current.step, bullets: result.summary?.bullets || [] })
    }
  } finally {
    sc.recipeProgress = null
  }
  return {
    ok: true,
    dataset: id,
    recipe: wanted,
    title: recipe.title || wanted,
    steps: done,
    completed: done.length,
    total: steps.length,
    suggestions: await suggestionsFor(runtime, id),
  }
}

export async function opExport(runtime, { dataset, what, path: outPath, actor } = {}) {
  const id = resolveDatasetId(runtime, dataset)
  const kind = String(what || 'h5ad')
  const allowed = ['h5ad', 'markers_csv', 'obs_csv', 'summary_json', 'protocol', 'methods', 'report', 'script']
  if (!allowed.includes(kind)) {
    throw Object.assign(new Error(`未知导出格式 "${what}"，可用：${allowed.join('/')}`), { code: 'bad_request' })
  }
  const engine = await ensureEngine(runtime)
  const started = Date.now()
  const result = await engine.request(
    'export',
    { dataset: id, what: kind, path: typeof outPath === 'string' && outPath.trim() ? outPath.trim() : undefined },
    { timeout: 180_000 },
  )
  logActivity(runtime, {
    actor: actorOf(actor), action: `export:${kind}`, dataset: id, ms: Date.now() - started,
    detail: `${result.file}（${Math.round((result.bytes || 0) / 1024)} KB）`,
  })
  return { ok: true, dataset: id, what: kind, file: result.file, bytes: result.bytes }
}

export async function opReset(runtime, { actor } = {}) {
  const sc = stateOf(runtime)
  const engine = runtime.engine
  if (engine && engine.running) {
    try {
      const result = await engine.request('list', {}, { timeout: 10_000 })
      for (const ds of result?.datasets || []) {
        try { await engine.request('close', { dataset: ds.id }, { timeout: 20_000 }) } catch { /* ignore */ }
      }
    } catch { /* ignore */ }
  }
  sc.datasets.clear()
  sc.figures.clear()
  sc.active = null
  sc.suggestCache = { key: '', value: [], at: 0 }
  sc.activity.length = 0
  logActivity(runtime, { actor: actorOf(actor), action: 'reset', detail: '清空工作台' })
  return { ok: true, active: null, datasets: [] }
}

export function opFigure(runtime, figureId) {
  const sc = stateOf(runtime)
  const id = path.basename(String(figureId || ''))
  if (!/^[A-Za-z0-9._-]+\.png$/.test(id)) throw Object.assign(new Error('非法的图片 id'), { code: 'bad_request' })
  const entry = sc.figures.get(id)
  if (!entry) throw Object.assign(new Error('图片不存在或已过期，请重新生成'), { code: 'not_found' })
  return entry
}

// ---------------------------------------------------------------------------
// Agent 工具
// ---------------------------------------------------------------------------

function renderStateForAgent(state) {
  const engine = state.engine || {}
  const lines = []
  const engineText = engine.ready
    ? `${engine.mode === 'scanpy' ? 'scanpy' : '内置引擎'} (Python ${engine.version || '?'})`
    : `未就绪${engine.error ? '：' + engine.error : ''}`
  lines.push(`单细胞工作台 · 引擎：${engineText}`)
  if (!state.datasets || !state.datasets.length) {
    lines.push('当前没有数据集。可 action=demo 载入内置示例数据，或 action=open 打开 .h5ad。')
  } else {
    lines.push(`数据集（${state.datasets.length} 个，active=${state.active}）：`)
    for (const d of state.datasets) {
      const done = (d.steps || []).length ? d.steps.join('→') : '无'
      lines.push(`- ${d.id}「${d.name}」${d.n_obs} 细胞 × ${d.n_vars} 基因；已完成：${done}`)
    }
  }
  if (state.suggestions && state.suggestions.length) {
    lines.push('建议下一步：' + state.suggestions.map((s) => s.title).join('；'))
  }
  return lines.join('\n')
}

function renderMarkersForAgent(groups, topN = 3) {
  const lines = []
  for (const g of groups) {
    const genes = (g.genes || []).slice(0, topN).map((x) => {
      const padj = x.pval_adj == null ? 'NA' : Number(x.pval_adj) < 1e-300 ? '<1e-300' : Number(x.pval_adj).toExponential(1)
      return `${x.gene}(logFC ${Number(x.logfc).toFixed(2)}, adj.P ${padj})`
    })
    lines.push(`簇 ${g.label ?? g.id}（n=${g.n_cells}）：${genes.join('、') || '无显著基因'}`)
  }
  return lines.join('\n')
}

function renderEmbeddingForAgent(result) {
  const cats = Array.isArray(result.categories) ? result.categories : []
  const head = `${result.layout} 坐标已计算：${result.n_points} 个点${result.downsampled ? '（已降采样）' : ''}，着色=${result.legend_title || ''}`
  if (!cats.length) return `${head}；数值范围 ${result.values ? `${Math.min(...result.values).toFixed(2)}–${Math.max(...result.values).toFixed(2)}` : 'n/a'}`
  const parts = cats.slice(0, 20).map((c) => `${c.label ?? c.id}:${c.n ?? '?'}`)
  return `${head}\n分组细胞数：${parts.join('，')}`
}

/**
 * 工具入口：参数与 UI 操作一一对应。
 * 结果保持紧凑（不回传上万点坐标），message 是给模型读的中文摘要。
 */
export async function runTool(runtime, args = {}) {
  const action = String(args.action || 'status')
  const actor = 'agent'
  const dataset = args.dataset ? String(args.dataset) : undefined
  switch (action) {
    case 'status': {
      const state = await opState(runtime)
      return {
        message: renderStateForAgent(state),
        active: state.active,
        datasets: state.datasets,
        suggestions: state.suggestions,
        engine: state.engine,
      }
    }
    case 'demo': {
      const result = await opDemo(runtime, { actor })
      const suggestions = await suggestionsFor(runtime, result.dataset.id)
      return {
        message: `已载入示例数据「${result.dataset.name}」：${result.dataset.n_obs} 细胞 × ${result.dataset.n_vars} 基因（id=${result.dataset.id}）。` +
          (suggestions.length ? `\n建议下一步：${suggestions.map((s) => s.title).join('；')}` : ''),
        dataset: result.dataset,
        suggestions,
      }
    }
    case 'open': {
      const result = await opOpen(runtime, {
        path: args.path, url: args.url, name: args.name, format: args.data_format || args.format,
        member: args.member, object: args.object, transpose: args.transpose,
        paths: args.paths, mergeDir: args.merge_dir, samples: args.sample_table, actor,
      })
      const suggestions = await suggestionsFor(runtime, result.dataset.id)
      return {
        message: `已导入「${result.dataset.name}」（格式 ${result.dataset.source_format || 'auto'}）：` +
          `${result.dataset.n_obs} 细胞 × ${result.dataset.n_vars} 基因（id=${result.dataset.id}）。` +
          (result.dataset.notes && result.dataset.notes.length ? `\n说明：${result.dataset.notes.join('；')}` : '') +
          (suggestions.length ? `\n建议下一步：${suggestions.map((s) => s.title).join('；')}` : ''),
        dataset: result.dataset,
        suggestions,
      }
    }
    case 'samples': {
      const op = String(args.op || 'list')
      const dataset = resolveDatasetId(runtime, args.dataset)
      let result
      if (op === 'list') result = await opObs(runtime, { dataset })
      else if (op === 'update') result = await opSamplesUpdate(runtime, { dataset, samples: args.samples, groups: args.groups, actor: 'agent' })
      else if (op === 'auto_group') result = await opSamplesAutoGroup(runtime, { dataset, pattern: args.group_regex ?? args.regex ?? args.pattern })
      else if (op === 'import_map') result = await opSamplesImportMap(runtime, { dataset, path: args.path })
      else if (op === 'join_obs') result = await opObsJoin(runtime, { dataset, path: args.path, key: args.key })
      else if (op === 'derive') result = await opObsDerive(runtime, { dataset, name: args.samples_op_name || args.name, source: args.source || 'cell_name', pattern: args.group_regex ?? args.regex ?? args.pattern })
      else if (op === 'set_group') result = await opObsSetGroup(runtime, { dataset, column: args.column })
      else if (op === 'remove') result = await opSamplesRemove(runtime, { dataset, samples: args.samples, actor: 'agent' })
      else throw Object.assign(new Error(`未知 samples op "${op}"，可用 list/update/auto_group/import_map/join_obs/derive/set_group/remove`), { code: 'bad_request' })
      if (op === 'remove') {
        return {
          message: `已移除样本 ${(result.removed_samples || []).join('、')}：共 ${result.removed_cells} 个细胞；聚类/嵌入等下游结果已重置，请重新运行降维聚类。`,
          dataset, removed_cells: result.removed_cells, removed_samples: result.removed_samples,
        }
      }
      const samples = (result.samples || []).map((item) => `- ${item.id} → 名称「${item.name}」大组「${item.group || '未设置'}」(${item.n_cells} 细胞)`).join('\n')
      return {
        message: `样本/分组（${op}）：${result.group_column ? '当前分组列=' + result.group_column : '尚未设置分组列'}` +
          (samples ? `\n样本表：\n${samples}` : '') +
          (result.matched !== undefined ? `\n匹配 ${result.matched}/${result.total}（${Math.round((result.match_rate || 0) * 100)}%）` : '') +
          (result.mapping ? `\n建议分组：${JSON.stringify(result.mapping)}` : ''),
        dataset, op, samples: result.samples || [], groups: result.groups || [],
        group_column: result.group_column || null, mapping: result.mapping || null,
        matched: result.matched, total: result.total,
      }
    }
    case 'scan': {
      const result = await opScan(runtime, { path: args.path, recursive: args.recursive, actor })
      const samples = Array.isArray(result.samples) ? result.samples : []
      const lines = samples.slice(0, 50).map((item) => `- ${item.name}${item.raw_name !== item.name ? `（原名 ${item.raw_name}）` : ''} · ${item.format_label} · ${item.detail || ''}`)
      return {
        message: `扫描完成：${result.root} 下识别到 ${samples.length} 个样本` +
          (result.errors && result.errors.length ? `；提示：${result.errors.join('；')}` : '') +
          (lines.length ? `\n${lines.join('\n')}` : '') +
          '\n下一步：把样本表（可改 name / 加 group）传给 action=open 的 sample_table 合并导入。',
        ...result,
      }
    }
    case 'geo': {
      const result = await opGeoLookup(runtime, { gse: args.gse || args.path })
      const samples = Array.isArray(result.samples) ? result.samples : []
      const lines = samples.slice(0, 40).map((item) => `- ${item.gsm}：${item.title || '（无标题）'}`)
      return {
        message: `GEO ${result.gse}：${result.title || ''}（${samples.length} 个样本）` +
          (lines.length ? `\n${lines.join('\n')}` : ''),
        ...result,
      }
    }
    case 'subset': {
      const result = await opSubset(runtime, { dataset: args.dataset, clusters: args.clusters, column: args.cluster_column, name: args.name, actor })
      return {
        message: `已提取小亚群：${result.dataset ? result.dataset.name : ''}（${result.subset_cells} 个细胞，簇 ${(result.kept_clusters || []).join('、')}）；已切换为当前数据集，可继续跑降维聚类。`,
        dataset: result.dataset, parent: result.parent, subset_cells: result.subset_cells,
      }
    }
    case 'inspect': {
      const result = await opInspect(runtime, {
        path: args.path, url: args.url, format: args.data_format || args.format, member: args.member, object: args.object,
      })
      const candidates = Array.isArray(result.candidates) ? result.candidates : []
      const lines = candidates.map((c) => `- ${c.member || c.object || c.name}${c.class ? '（' + c.class + '）' : ''}${c.n_obs ? ` ${c.n_obs}×${c.n_vars}` : ''}`)
      return {
        message: `格式探测：${result.format || '未知'}` +
          (result.n_obs ? `，约 ${result.n_obs} 细胞 × ${result.n_vars} 基因` : '') +
          (result.notes && result.notes.length ? `\n说明：${result.notes.join('；')}` : '') +
          (lines.length ? `\n候选对象/文件：\n${lines.join('\n')}\n（用 open 时带 object 或 member 指定）` : ''),
        ...result,
      }
    }
    case 'analyses': {
      const result = await opAnalyses(runtime, { dataset })
      const native = result.analyses.filter((item) => item.native)
      const usable = native.filter((item) => item.available)
      const external = result.analyses.filter((item) => !item.native)
      const lines = native.map((item) => `- ${item.id}（${item.group}）${item.available ? '可用' : '不可用：' + (item.reason || '')}`)
      return {
        message: `高级分析：${usable.length}/${native.length} 项当前可运行，另有 ${external.length} 项需要外部工具（不作为原生分析运行）。\n` +
          lines.join('\n') + '\n用 action=analysis + analysis=<id> 运行。',
        groups: result.groups,
        dataset: result.dataset,
        analyses: result.analyses.map((item) => ({
          id: item.id, title: item.title, group: item.group, method: item.method || '',
          native: Boolean(item.native), available: Boolean(item.available),
          reason: item.reason || null, params: Array.isArray(item.params) ? item.params : [],
        })),
        diagnostics: result.diagnostics,
      }
    }
    case 'analysis': {
      const result = await opAnalysis(runtime, { dataset, id: args.analysis || args.id, params: args.params, actor })
      const bullets = (result.summary?.bullets || []).map((line) => `- ${line}`)
      const tables = (result.tables || []).map((table) => `表「${table.title || table.id}」${(table.rows || []).length} 行`)
      const figures = (result.figures || []).map((figure) => figure.file).filter(Boolean)
      return {
        message: `已运行「${result.title}」（${result.method || ''}）` +
          (bullets.length ? '\n' + bullets.join('\n') : '') +
          (tables.length ? '\n' + tables.join('；') : '') +
          (result.obs_merged && result.obs_merged.length ? `\n已写回 obs 列：${result.obs_merged.join('、')}` : '') +
          (figures.length ? `\n图文件（可用 read_image 查看）：\n${figures.join('\n')}` : '') +
          ((result.warnings || []).length ? `\n注意：${(result.warnings || []).join('；')}` : ''),
        id: result.id,
        title: result.title,
        method: result.method,
        summary: result.summary,
        tables: (result.tables || []).map((table) => ({
          id: table.id, title: table.title, columns: table.columns, rows: (table.rows || []).slice(0, 20),
        })),
        figures: result.figures,
        obs_merged: result.obs_merged,
        warnings: result.warnings,
      }
    }
    case 'preview': {
      const result = await opPreview(runtime, { dataset })
      return {
        message: `${result.dataset}：${result.n_obs} 细胞 × ${result.n_vars} 基因，已完成 ${(result.steps || []).length} 步。\n` +
          `部分基因：${(result.genes || []).slice(0, 20).join(', ')}\n` +
          `obs 列：${(result.obs_columns || []).slice(0, 12).map((c) => c.key).join(', ')}`,
        genes: (result.genes || []).slice(0, 50),
        obs_columns: result.obs_columns,
        steps: result.steps,
      }
    }
    case 'run': {
      const result = await opRun(runtime, { dataset, step: args.step, params: args.params, actor })
      const summary = result.summary || {}
      const bullets = (summary.bullets || []).map((b) => `- ${b}`).join('\n')
      const next = (result.suggestions || []).map((s) => s.title).join('；')
      return {
        message: `已执行「${summary.title || result.step}」${bullets ? '\n' + bullets : ''}` +
          (next ? `\n建议下一步：${next}` : ''),
        step: result.step,
        summary: summary.metrics || {},
        steps: result.steps,
        suggestions: result.suggestions,
      }
    }
    case 'embedding': {
      const result = await opEmbedding(runtime, { dataset, color_by: args.color_by })
      return {
        message: renderEmbeddingForAgent(result) + '\n（坐标只在工作台画布渲染，不回传原始点；需要看图请 action=plot kind=embedding）',
        layout: result.layout,
        n_points: result.n_points,
        categories: result.categories,
      }
    }
    case 'markers': {
      const result = await opMarkers(runtime, { dataset, top_n: args.top_n, method: args.method, actor })
      return {
        message: `各簇 top 标记基因（${result.groups.length} 簇）：\n` + renderMarkersForAgent(result.groups, 5),
        groups: result.groups.map((g) => ({ id: g.id, label: g.label, n_cells: g.n_cells, genes: (g.genes || []).slice(0, 8) })),
      }
    }
    case 'audit': {
      const result = await opAudit(runtime, { dataset, actor })
      const mark = (level) => (level === 'ok' ? 'OK ' : level === 'warn' ? 'WARN' : 'FAIL')
      const lines = result.checks.map((c) => `[${mark(c.level)}] ${c.title}：${c.detail}${c.hint ? ' → ' + c.hint : ''}`)
      return {
        message: `自检得分 ${result.score}/100（${result.level}）\n${lines.join('\n')}`,
        score: result.score,
        level: result.level,
        checks: result.checks,
      }
    }
    case 'autoannotate': {
      const result = await opAutoannotate(runtime, { dataset, actor })
      if (!result.suggestions.length) {
        return { message: '暂时无法给出自动注释建议（需要先完成聚类；参考面板基因可能不在该数据集中）。', suggestions: [] }
      }
      const lines = result.suggestions.map((s) => {
        const ev = (s.evidence || []).join('/')
        const alt = s.alternative ? `；备选 ${s.alternative}` : ''
        return `簇 ${s.cluster} → ${s.label}（score ${Number(s.score).toFixed(2)}，margin ${Number(s.margin ?? 0).toFixed(2)}；证据 ${ev}${alt}）`
      })
      return {
        message:
          `参考 marker 面板 ${result.panel} 的自动注释建议（只给建议，未写入 obs）：\n${lines.join('\n')}\n` +
          '确认后调用 action=run step=annotate params={"mapping":{"0":"T cell",...}} 写回细胞类型。',
        panel: result.panel,
        suggestions: result.suggestions,
      }
    }
    case 'resolution_scan': {
      const result = await opResolutionScan(runtime, { dataset, resolutions: args.resolutions, actor })
      const cur = result.current ? `当前 resolution=${result.current.resolution}（${result.current.n_clusters} 簇）` : '当前未聚类'
      const lines = result.results.map((r) => `resolution ${r.resolution} → ${r.n_clusters} 簇，ARI vs 当前 ${Number(r.ari_vs_current).toFixed(3)}`)
      return {
        message: `聚类分辨率扫描（不改变当前聚类）：${cur}\n${lines.join('\n')}\n结论建议：选 ARI 明显下降前、簇数符合生物学预期的分辨率。`,
        current: result.current,
        results: result.results,
      }
    }
    case 'recipe': {
      const result = await opRecipe(runtime, { dataset, id: args.recipe || args.id || 'standard', actor })
      const lines = result.steps.map((s, i) => `${i + 1}. ${s.title}${s.bullets.length ? '：' + s.bullets.join('；') : ''}`)
      const next = (result.suggestions || []).map((s) => s.title).join('；')
      return {
        message: `已运行预设流程「${result.title}」（${result.completed}/${result.total} 步完成）\n${lines.join('\n')}` +
          (next ? `\n建议下一步：${next}` : ''),
        recipe: result.recipe,
        steps: result.steps,
        suggestions: result.suggestions,
      }
    }
    case 'plot': {
      const result = await opPlot(runtime, { dataset, kind: args.kind, params: args.params, actor })
      return {
        message: `已生成 ${result.figure.kind} 图（${result.figure.width}×${result.figure.height}px），文件：${result.figure.file}\n` +
          '可在工作台对应 Tab 查看；需要模型自己看图时用 read_image 打开该文件。',
        figure: result.figure,
      }
    }
    case 'export': {
      const result = await opExport(runtime, { dataset, what: args.what, path: args.out_path, actor })
      return { message: `已导出 ${result.what} → ${result.file}（${Math.round((result.bytes || 0) / 1024)} KB）`, ...result }
    }
    case 'close': {
      const result = await opClose(runtime, { dataset, actor })
      return { message: `已关闭数据集 ${dataset || ''}；剩余 ${result.datasets.length} 个。`, ...result }
    }
    case 'reset': {
      const result = await opReset(runtime, { actor })
      return { message: '工作台已清空。', ...result }
    }
    default:
      throw Object.assign(
        new Error(`未知 action "${action}"，可用：status/demo/open/preview/run/embedding/markers/plot/export/close/reset`),
        { code: 'bad_request' },
      )
  }
}

/** 工具输出必须是 lossless JSON：去掉 undefined、把 NaN/Infinity 转 null。 */
function sanitizeValue(value, depth = 0) {
  if (depth > 16) return null
  if (value === undefined) return undefined
  if (typeof value === 'number') return Number.isFinite(value) ? value : null
  if (typeof value === 'bigint') return Number(value)
  if (typeof value === 'function' || typeof value === 'symbol') return null
  if (Array.isArray(value)) {
    return value.map((item) => {
      const next = sanitizeValue(item, depth + 1)
      return next === undefined ? null : next
    })
  }
  if (value && typeof value === 'object') {
    const out = {}
    for (const [key, item] of Object.entries(value)) {
      const next = sanitizeValue(item, depth + 1)
      if (next !== undefined) out[key] = next
    }
    return out
  }
  return value
}

function buildToolDefinition(hostRuntime) {
  const def = {
    name: 'singlecell',
    description:
      '单细胞分析工作台（本机已安装的 Host 插件，带可视化页面与计算引擎）。用它打开或分析单细胞数据（.h5ad）、跑标准流程、取聚类/标记基因/嵌入图结果。' +
      'action: status=查看工作台状态（已载入数据集、已跑步骤、建议下一步）; ' +
      'demo=载入内置示例数据（无需文件）; ' +
      'open=导入数据（path 本地路径 或 url 直链；支持 h5ad、10x .h5、10x mtx 三件套/目录、loom、csv/tsv/txt 表达矩阵、GEO series matrix、.rds（Seurat/SCE）、.RData/.rda、zip/tar.gz 压缩包；可选 data_format/member/object/transpose）; ' +
      'inspect=只探测格式与维度（大文件/压缩包有歧义时先用它）; ' +
      'preview=查看基因名与 obs 列; run=运行一个流程步骤（step 必填，可选 params）; ' +
      'embedding=计算二维嵌入坐标并返回分组统计; markers=返回各簇 top 标记基因; ' +
      'autoannotate=用内置参考 marker 面板给出每簇细胞类型建议（只给建议，采用需 run step=annotate）; ' +
      'audit=确定性自检（流程完整性、线粒体比例、簇数、双细胞风险等，返回 score 与逐项问题）; ' +
      'resolution_scan=扫描多个聚类分辨率并返回簇数与 ARI（参数对比，不改当前聚类）; ' +
      'recipe=一键跑预设流程（standard/quick/qc）; ' +
      'samples=样本/分组管理（op=list 看样本表; op=update 改样本名或设大组; op=auto_group 按 regex 正则批量分组; op=import_map 读样本→大组对照 CSV; op=join_obs 把细胞级 metadata CSV join 进 obs; op=derive 用 regex 从细胞名派生列; op=set_group 设为分组列）——GEO 每 GSM 一个文件时先用 paths 合并导入再设分组; ' +
      'analyses=列出全部高级分析（基本分析/组间对比/可视化/通路/基因集/轨迹/异质性/通讯/转录调控/虚拟敲除）及可用性; ' +
      'analysis=运行一个高级分析（analysis 传 id，params 传参数；结果含表/图，部分分析会写回 obs 列）; ' +
      'plot=出图（kind=qc_violin/embedding/dotplot/heatmap/cluster_bar/pca_variance/qc_scatter，返回 PNG 绝对路径）; ' +
      'export=导出（what=h5ad/markers_csv/obs_csv/summary_json/protocol/methods/report/script，可用 out_path 指定路径）; ' +
      'close=关闭数据集; reset=清空工作台。' +
      '流程步骤顺序：qc→filter→normalize→hvg→pca→neighbors→cluster→embed→markers→annotate（每一步依赖前一步；cluster 之后才能 markers/annotate）。' +
      '用户提到单细胞、聚类、UMAP/t-SNE、标记基因、细胞类型注释、h5ad 时优先用本工具驱动工作台（与页面按钮同一套实现），不要自己写 Python 脚本重跑同一套分析。' +
      '给结论前建议先 action=audit，注释前先 action=autoannotate 看参考面板建议。不确定状态时先调 action=status。',
    parameters: {
      action: {
        type: 'string',
        enum: ['status', 'demo', 'scan', 'geo', 'open', 'inspect', 'preview', 'run', 'embedding', 'markers', 'autoannotate', 'audit', 'resolution_scan', 'recipe', 'analyses', 'analysis', 'samples', 'subset', 'plot', 'export', 'close', 'reset'],
        description: '要执行的操作',
      },
      dataset: { type: 'string', description: '数据集 id（如 ds1）；省略则用当前活动数据集' },
      path: { type: 'string', description: '本地文件/目录绝对路径（action=open/inspect/scan 用；scan 传 GEO 解压后的数据根目录）' },
      url: { type: 'string', description: 'http(s)/ftp 直链（GEO FTP 链接等；action=open/inspect 用）' },
      // 注意：参数名不能用 `format`/`pattern` 这类 JSON Schema 保留关键字（要求 string 值）。
      // Host 把 parameters 原样投给模型侧，撞上保留字会被 API 判为非法 schema（两次教训）。
      data_format: {
        type: 'string',
        description: '强制格式（可选，默认自动）：h5ad|10x_h5|10x_mtx|loom|matrix_text|geo_series|rds|rdata|archive',
      },
      member: { type: 'string', description: '压缩包内的相对路径（archive 有多个候选时用）' },
      object: { type: 'string', description: 'RData 内的对象名（多个对象候选时用）' },
      transpose: { type: 'boolean', description: '表达矩阵方向（true=基因×细胞 转置成 细胞×基因）' },
      analysis: { type: 'string', description: 'action=analysis 的分析 id（如 cell_cycle/de_group/pathway_ora/integration/trajectory_dpt/lr_communication）' },
      paths: { type: 'array', description: 'action=open 的多文件合并导入：本地文件路径数组（每个文件=一个样本，自动加 sample 列，基因取交集）', items: { type: 'string' } },
      merge_dir: { type: 'boolean', description: 'action=open 时把 path 指向的目录下每个数据文件当作一个样本合并导入' },
      recursive: { type: 'boolean', description: 'action=scan 是否递归扫描子目录（默认只扫一层）' },
      gse: { type: 'string', description: 'action=geo 的 GSE 编号（如 GSE246662），返回每个 GSM 的说明辅助分组' },
      sample_table: { type: 'array', description: 'action=open 的样本表（scan 结果勾选/改名/分组后）：[{path, name, group, trio_prefix, format}]', items: { type: 'object' } },
      clusters: { type: 'array', description: 'action=subset 要提取的簇号列表（如 ["3","7"]）', items: { type: 'string' } },
      cluster_column: { type: 'string', description: 'action=subset 依据的 obs 列名（默认 clusters）' },
      op: { type: 'string', description: 'action=samples 的子操作：list | update | auto_group | import_map | join_obs | derive | set_group | remove（按样本真正剔除细胞）' },
      samples: { type: 'array', description: 'action=samples op=update 的样本列表：[{id, name, group}]（id=源文件名，name=自定义样本名，group=大组）', items: { type: 'object' } },
      groups: { type: 'array', description: 'action=samples op=update 的大组列表（可选，如 ["young","old"]）', items: { type: 'string' } },
      column: { type: 'string', description: 'action=samples op=set_group 要设为分组列的 obs 列名' },
      group_regex: { type: 'string', description: 'action=samples op=auto_group 的文件名/样本名正则（第 1 个捕获组=大组）；op=derive 时用于细胞名正则' },
      samples_op_name: { type: 'string', description: 'action=samples op=derive 的新列名' },
      name: { type: 'string', description: '数据集显示名（可选）' },
      step: {
        type: 'string',
        enum: STEP_IDS,
        description: '流程步骤（action=run 必填）',
      },
      params: {
        type: 'object',
        description:
          '步骤参数。qc/filter: min_genes,max_genes,max_pct_mt,min_counts,min_cells；normalize: target_sum；hvg: n_top_genes；' +
          'pca: n_pcs；neighbors: n_neighbors,metric；cluster: resolution,method；embed: layout(tsne|umap|pca),perplexity；' +
          'markers: top_n,method；annotate: mapping({"0":"T cell"})',
        additionalProperties: true,
      },
      color_by: {
        type: 'object',
        description: 'action=embedding 的着色：{"kind":"cluster"|"obs"|"gene","key":"clusters"}',
        additionalProperties: true,
      },
      kind: { type: 'string', enum: PLOT_KINDS, description: 'action=plot 的图类型' },
      top_n: { type: 'number', description: '每簇返回的标记基因数（默认 10）' },
      method: { type: 'string', description: 'markers 统计方法（默认 t-test）' },
      what: {
        type: 'string',
        enum: ['h5ad', 'markers_csv', 'obs_csv', 'summary_json', 'protocol', 'methods', 'report', 'script'],
        description: 'action=export 的导出格式',
      },
      out_path: { type: 'string', description: 'action=export 的目标文件绝对路径（可选）' },
      recipe: { type: 'string', description: 'action=recipe 的预设 id：standard（推荐全流程）/quick/qc' },
      resolutions: { type: 'string', description: 'action=resolution_scan 的分辨率，逗号分隔，如 "0.4,0.8,1.2"' },
    },
    output: {
      schema: { type: 'object' },
      render: (_args, value) => {
        const text = value && typeof value.message === 'string' ? value.message : JSON.stringify(value)
        return [{ type: 'text', text: String(text) }]
      },
    },
    execute: async (args) => {
      const fresh = await import(new URL(`./impl.mjs?rev=${Date.now()}-${Math.random().toString(36).slice(2)}`, import.meta.url).href)
      return sanitizeValue(await fresh.runTool(hostRuntime, args || {}))
    },
  }
  // 防呆：参数名撞 JSON Schema 保留关键字时，API 会把整个 tools 数组判非法（format、pattern 各炸过一次）。
  // 注册期就抛错并写日志：只丢本工具，不拖垮宿主里所有 agent 请求。
  const reservedKeywords = new Set(['type', 'enum', 'const', 'multipleOf', 'maximum', 'exclusiveMaximum', 'minimum', 'exclusiveMinimum', 'maxLength', 'minLength', 'pattern', 'maxItems', 'minItems', 'uniqueItems', 'maxContains', 'minContains', 'maxProperties', 'minProperties', 'required', 'dependentRequired', 'properties', 'patternProperties', 'additionalProperties', 'propertyNames', 'items', 'prefixItems', 'additionalItems', 'contains', 'title', 'description', 'default', 'deprecated', 'readOnly', 'writeOnly', 'examples', 'format', 'contentEncoding', 'contentMediaType', 'contentSchema', 'if', 'then', 'else', 'allOf', 'anyOf', 'oneOf', 'not'])
  const clash = Object.keys(def.parameters).filter((key) => reservedKeywords.has(key))
  if (clash.length) throw new Error(`[dsh-singlecell-workbench] 参数名 ${clash.join('、')} 是 JSON Schema 保留关键字，请改名并 bump TOOL_SCHEMA_VERSION（format、pattern 都栽过这个坑）`)
  return def
}

// 工具 schema 只在注册时快照一次。改了参数定义必须 bump 这个版本号，
// 让下一次 /singlecell/* 请求先注销旧定义再注册新定义（无需重启宿主）。
const TOOL_SCHEMA_VERSION = 4 // v4: regex -> group_regex（参数名避开 JSON Schema 保留字；v3 regex、v2 format -> data_format）

export function ensureToolRegistered(runtime) {
  if (runtime.toolRegistered && runtime.toolSchemaVersion === TOOL_SCHEMA_VERSION) return runtime.toolRegistered
  try {
    const ctx = runtime.ctx
    const tools = ctx && typeof ctx.get === 'function' ? ctx.get('tools') : ctx && ctx.tools
    if (!tools || typeof tools.register !== 'function') {
      runtime.toolRegistered = false
      return false
    }
    if (runtime.toolRegistered) {
      // 拿不到 disposer 就别动注册表，否则同名注册会抛错、反而把工具弄丢。
      if (typeof runtime.toolRegistered !== 'function') return runtime.toolRegistered
      try { runtime.toolRegistered() } catch (error) {
        try { console.warn('[dsh-singlecell-workbench] old tool dispose failed:', error && error.message ? error.message : error) } catch { /* ignore */ }
      }
      runtime.toolRegistered = false
    }
    const dispose = tools.register(buildToolDefinition(runtime))
    runtime.toolRegistered = typeof dispose === 'function' ? dispose : true
    runtime.toolSchemaVersion = TOOL_SCHEMA_VERSION
    return true
  } catch (error) {
    runtime.toolRegistered = false
    try { console.warn('[dsh-singlecell-workbench] agent tool unavailable:', error && error.message ? error.message : error) } catch { /* ignore */ }
    return false
  }
}

// ---------------------------------------------------------------------------
// Prompt section：让 Agent 知道工作台存在、当前在分析什么
// ---------------------------------------------------------------------------

export function renderPromptSection(runtime) {
  const sc = runtime.sc
  if (!sc || !sc.datasets || sc.datasets.size === 0) return ''
  const active = sc.active ? sc.datasets.get(sc.active) : null
  const lines = ['单细胞分析工作台（Host 插件 dsh-singlecell-workbench）已在使用中。']
  if (active) {
    const done = (active.steps || []).length ? active.steps.join('、') : '无'
    lines.push(`当前数据集：${active.name || active.id}（${active.n_obs} 细胞 × ${active.n_vars} 基因），已完成步骤：${done}。`)
  }
  lines.push(
    '用户提到单细胞/聚类/UMAP/标记基因/细胞类型时，优先调用 singlecell 工具操作这个工作台（与页面按钮同一套实现），' +
      '不要另写脚本重跑同一套分析；图可用 read_image 查看返回的 PNG 路径。',
  )
  return lines.join('\n')
}

export function ensurePromptSection(runtime) {
  if (runtime.promptRegistered) return runtime.promptRegistered
  try {
    const ctx = runtime.ctx
    const systemPrompt = ctx && typeof ctx.get === 'function' ? ctx.get('systemPrompt') : ctx && ctx.systemPrompt
    if (!systemPrompt || typeof systemPrompt.section !== 'function') {
      runtime.promptRegistered = false
      return false
    }
    const dispose = systemPrompt.section({
      name: 'dsh-singlecell-workbench',
      order: 760,
      text: () => renderPromptSection(runtime),
    })
    runtime.promptRegistered = typeof dispose === 'function' ? dispose : true
    return true
  } catch (error) {
    runtime.promptRegistered = false
    try { console.warn('[dsh-singlecell-workbench] prompt section unavailable:', error && error.message ? error.message : error) } catch { /* ignore */ }
    return false
  }
}

// ---------------------------------------------------------------------------
// 路由分发
// ---------------------------------------------------------------------------

async function handleUpload(runtime, req, res) {
  const started = Date.now()
  const filename = sanitizeName(req.headers['x-dsh-filename'] || 'upload.h5ad')
  const engine = await ensureEngine(runtime)
  const dest = path.join(engine.dirs.uploads, `${Date.now().toString(36)}-${randomUUID().slice(0, 8)}-${filename}`)
  let size = 0
  try {
    await pipeline(
      req,
      async function* count(source) {
        for await (const chunk of source) {
          size += chunk.length
          if (size > MAX_UPLOAD_BYTES) throw Object.assign(new Error('文件超过 8 GiB 上限'), { code: 'bad_request' })
          yield chunk
        }
      },
      createWriteStream(dest),
    )
  } catch (error) {
    try { unlinkSync(dest) } catch { /* ignore */ }
    logActivity(runtime, { actor: 'user', action: 'upload', detail: `${filename} → ${error.message}`, ok: false, ms: Date.now() - started })
    throw error
  }
  return opOpen(runtime, { path: dest, name: filename, actor: 'user' })
}

async function serveFigure(runtime, url, res) {
  const id = url.searchParams.get('id')
  const entry = opFigure(runtime, id)
  let stat
  try { stat = statSync(entry.file) } catch {
    throw Object.assign(new Error('图片文件已被清理，请重新生成'), { code: 'not_found' })
  }
  res.writeHead(200, {
    'content-type': 'image/png',
    'cache-control': 'private, max-age=3600',
    'content-length': String(stat.size),
  })
  await pipeline(createReadStream(entry.file), res)
}

export async function handle(runtime, req, res) {
  ensureToolRegistered(runtime)
  ensurePromptSection(runtime)
  if (!hostIsLoopback(req)) {
    return sendJson(res, 403, { ok: false, code: 'forbidden', error: '拒绝：只接受本机（回环地址）请求' })
  }
  const url = new URL(req.url || '/', 'http://127.0.0.1')
  if (!url.pathname.startsWith(PREFIX)) {
    return sendJson(res, 404, { ok: false, code: 'not_found', error: '未知路由' })
  }
  const route = url.pathname.slice(PREFIX.length) || '/'
  const method = (req.method || 'GET').toUpperCase()
  const key = `${method} ${route}`
  try {
    switch (key) {
      case 'GET /health':
        return sendJson(res, 200, await opHealth(runtime))
      case 'GET /state':
        return sendJson(res, 200, await opState(runtime))
      case 'GET /activity':
        return sendJson(res, 200, { ok: true, entries: activityOf(runtime, Number(url.searchParams.get('limit')) || 80) })
      case 'GET /figure':
        return await serveFigure(runtime, url, res)
      case 'POST /dataset/open': {
        const body = await readJson(req)
        return sendJson(res, 200, await opOpen(runtime, {
          path: body.path, url: body.url, name: body.name, format: body.format,
          member: body.member, object: body.object, transpose: body.transpose,
          paths: body.paths, mergeDir: body.merge_dir, samples: body.samples,
          actor: actorOf(body.actor),
        }))
      }
      case 'GET /formats':
        return sendJson(res, 200, await opFormats(runtime))
      case 'POST /pick': {
        const body = await readJson(req)
        return sendJson(res, 200, await opPickPath(runtime, { kind: body.kind, prompt: body.prompt, dryRun: body.dry_run }))
      }
      case 'POST /dataset/subset': {
        const body = await readJson(req)
        return sendJson(res, 200, await opSubset(runtime, { dataset: body.dataset, clusters: body.clusters, column: body.column || body.cluster_column, name: body.name, actor: actorOf(body.actor) }))
      }
      case 'POST /scan': {
        const body = await readJson(req)
        return sendJson(res, 200, await opScan(runtime, { path: body.path || body.root, recursive: body.recursive, actor: actorOf(body.actor) }))
      }
      case 'POST /geo': {
        const body = await readJson(req)
        return sendJson(res, 200, await opGeoLookup(runtime, { gse: body.gse || body.id }))
      }
      case 'POST /import/inspect': {
        const body = await readJson(req)
        return sendJson(res, 200, await opInspect(runtime, {
          path: body.path, url: body.url, format: body.format, member: body.member, object: body.object,
        }))
      }
      case 'GET /obs': {
        const id = url.searchParams.get('dataset') || undefined
        return sendJson(res, 200, await opObs(runtime, { dataset: id }))
      }
      case 'POST /samples/update': {
        const body = await readJson(req)
        return sendJson(res, 200, await opSamplesUpdate(runtime, {
          dataset: body.dataset, samples: body.samples, groups: body.groups,
          sampleColumn: body.sample_column, groupColumn: body.group_column, actor: actorOf(body.actor),
        }))
      }
      case 'POST /samples/remove': {
        const body = await readJson(req)
        return sendJson(res, 200, await opSamplesRemove(runtime, { dataset: body.dataset, samples: body.samples, actor: actorOf(body.actor) }))
      }
      case 'POST /samples/auto-group': {
        const body = await readJson(req)
        return sendJson(res, 200, await opSamplesAutoGroup(runtime, { dataset: body.dataset, pattern: body.pattern }))
      }
      case 'POST /samples/import-map': {
        const body = await readJson(req)
        return sendJson(res, 200, await opSamplesImportMap(runtime, {
          dataset: body.dataset, path: body.path, sampleColumn: body.sample_column, groupColumn: body.group_column,
        }))
      }
      case 'POST /obs/join': {
        const body = await readJson(req)
        return sendJson(res, 200, await opObsJoin(runtime, { dataset: body.dataset, path: body.path, key: body.key }))
      }
      case 'POST /obs/derive': {
        const body = await readJson(req)
        return sendJson(res, 200, await opObsDerive(runtime, {
          dataset: body.dataset, name: body.name, source: body.source, pattern: body.pattern, group: body.group, value: body.value,
        }))
      }
      case 'POST /obs/set-group': {
        const body = await readJson(req)
        return sendJson(res, 200, await opObsSetGroup(runtime, { dataset: body.dataset, column: body.column }))
      }
      case 'GET /external-tools':
        return sendJson(res, 200, { ok: true, ...(await (async () => {
          const engine = await ensureEngine(runtime)
          const sc = stateOf(runtime)
          return engine.request('external_tools', sc.active ? { dataset: sc.active } : {}, { timeout: 300_000 })
        })()) })
      case 'GET /analyses':
        return sendJson(res, 200, await opAnalyses(runtime, { dataset: url.searchParams.get('dataset') || undefined }))
      case 'POST /analysis': {
        const body = await readJson(req)
        return sendJson(res, 200, await opAnalysis(runtime, {
          dataset: body.dataset, id: body.id || body.analysis, params: body.params, actor: actorOf(body.actor),
        }))
      }
      case 'POST /dataset/demo': {
        const body = await readJson(req)
        return sendJson(res, 200, await opDemo(runtime, { actor: actorOf(body.actor) }))
      }
      case 'POST /dataset/upload':
        return sendJson(res, 200, await handleUpload(runtime, req, res))
      case 'POST /dataset/select': {
        const body = await readJson(req)
        return sendJson(res, 200, await opSelect(runtime, { dataset: body.dataset }))
      }
      case 'POST /dataset/close': {
        const body = await readJson(req)
        return sendJson(res, 200, await opClose(runtime, { dataset: body.dataset, actor: actorOf(body.actor) }))
      }
      case 'POST /preview': {
        const body = await readJson(req)
        return sendJson(res, 200, await opPreview(runtime, { dataset: body.dataset }))
      }
      case 'POST /run': {
        const body = await readJson(req)
        return sendJson(res, 200, await opRun(runtime, { dataset: body.dataset, step: body.step, params: body.params, actor: actorOf(body.actor) }))
      }
      case 'POST /embedding': {
        const body = await readJson(req)
        return sendJson(res, 200, await opEmbedding(runtime, { dataset: body.dataset, color_by: body.color_by }))
      }
      case 'POST /plot': {
        const body = await readJson(req)
        return sendJson(res, 200, await opPlot(runtime, { dataset: body.dataset, kind: body.kind, params: body.params, actor: actorOf(body.actor) }))
      }
      case 'POST /markers': {
        const body = await readJson(req)
        return sendJson(res, 200, await opMarkers(runtime, { dataset: body.dataset, top_n: body.top_n, method: body.method, actor: actorOf(body.actor) }))
      }
      case 'POST /audit': {
        const body = await readJson(req)
        return sendJson(res, 200, await opAudit(runtime, { dataset: body.dataset, actor: actorOf(body.actor) }))
      }
      case 'POST /autoannotate': {
        const body = await readJson(req)
        return sendJson(res, 200, await opAutoannotate(runtime, { dataset: body.dataset, actor: actorOf(body.actor) }))
      }
      case 'POST /resolution-scan': {
        const body = await readJson(req)
        return sendJson(res, 200, await opResolutionScan(runtime, { dataset: body.dataset, resolutions: body.resolutions, actor: actorOf(body.actor) }))
      }
      case 'POST /recipe': {
        const body = await readJson(req)
        return sendJson(res, 200, await opRecipe(runtime, { dataset: body.dataset, id: body.id, actor: actorOf(body.actor) }))
      }
      case 'POST /export': {
        const body = await readJson(req)
        return sendJson(res, 200, await opExport(runtime, { dataset: body.dataset, what: body.what, path: body.out_path || body.path, actor: actorOf(body.actor) }))
      }
      case 'POST /reset': {
        const body = await readJson(req)
        return sendJson(res, 200, await opReset(runtime, { actor: actorOf(body.actor) }))
      }
      default:
        return sendJson(res, 404, { ok: false, code: 'not_found', error: `未知路由 ${key}` })
    }
  } catch (error) {
    const failure = toFailure(error)
    if (failure.status >= 500) {
      try { console.error('[dsh-singlecell-workbench]', key, error && (error.stack || error.message || error)) } catch { /* ignore */ }
    }
    try {
      if (!res.headersSent) sendJson(res, failure.status, failure.body)
      else res.end()
    } catch { /* ignore */ }
  }
}

export { PREFIX, STEPS as STEP_DEFS }
