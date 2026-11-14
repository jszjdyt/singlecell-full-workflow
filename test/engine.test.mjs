/**
 * dsh-singlecell-workbench — 端到端联调测试（不装插件、不启 Cordis）
 *
 *   node test/engine.test.mjs
 *
 * 直接构造一个假 runtime，调 lib/impl.mjs 的 op* 函数（= UI 与 Agent 共用实现），
 * 走真实 Python worker：open demo → 全流程 → 嵌入坐标 → 出图 → 标记基因 → 导出 → 清理。
 *
 * 通过退出码报告：0 全过，1 有失败。
 */

import assert from 'node:assert/strict'
import { existsSync, mkdtempSync, readFileSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'

const here = new URL('.', import.meta.url)
const workdir = mkdtempSync(path.join(tmpdir(), 'dsh-singlecell-workbench-test-'))
const results = []
const t0 = Date.now()

function check(name, fn) {
  const started = Date.now()
  try {
    const value = fn()
    results.push({ name, ok: true, ms: Date.now() - started })
    console.log(`PASS  ${name}  (${Date.now() - started}ms)`)
    return value
  } catch (error) {
    results.push({ name, ok: false, ms: Date.now() - started, error })
    console.log(`FAIL  ${name}  (${Date.now() - started}ms)\n      ${error && error.message ? error.message : error}`)
    return undefined
  }
}

async function checkAsync(name, fn) {
  const started = Date.now()
  try {
    const value = await fn()
    results.push({ name, ok: true, ms: Date.now() - started })
    console.log(`PASS  ${name}  (${Date.now() - started}ms)`)
    return value
  } catch (error) {
    results.push({ name, ok: false, ms: Date.now() - started, error })
    console.log(`FAIL  ${name}  (${Date.now() - started}ms)\n      ${error && error.message ? error.message : error}`)
    return undefined
  }
}

const impl = await import(new URL('../lib/impl.mjs', import.meta.url).href)
const engineMod = await import(new URL('../lib/engine.mjs', import.meta.url).href)

console.log('workdir =', workdir)
console.log('engine script =', engineMod.showEngineScript ? engineMod.showEngineScript() : '(n/a)')

// 假 runtime：只有 config 与可选 ctx，足够驱动 op*
const runtime = { config: { workdir }, ctx: null }

let datasetId = null
let figureFile = null

await checkAsync('engine: 探测 Python + 启动 worker + ping 握手', async () => {
  const engine = await engineMod.createEngine({ workdir })
  const info = await engine.start()
  assert.ok(info && info.python, '应返回 python 路径')
  assert.ok(['builtin', 'scanpy'].includes(info.mode), `mode 非法：${info.mode}`)
  console.log(`      python=${info.python} version=${info.version} mode=${info.mode}`)
  assert.ok(info.packages && info.packages.numpy, 'packages.numpy 缺失')
  // 交给 runtime 复用同一个引擎实例
  runtime.engine = engine
})

await checkAsync('health: 引擎就绪状态', async () => {
  const health = await impl.opHealth(runtime)
  assert.equal(health.ok, true)
  assert.equal(health.engine.ready, true, `引擎未就绪：${health.engine.error || ''}`)
  assert.ok(health.engine.python, 'python 路径缺失')
})

await checkAsync('demo: 载入内置示例数据', async () => {
  const result = await impl.opDemo(runtime, { actor: 'user' })
  assert.equal(result.ok, true)
  const ds = result.dataset
  assert.ok(ds && ds.id, 'dataset.id 缺失')
  assert.ok(ds.n_obs >= 2000, `n_obs 太小：${ds.n_obs}`)
  assert.ok(ds.n_vars >= 1000, `n_vars 太小：${ds.n_vars}`)
  datasetId = ds.id
  console.log(`      dataset=${ds.id} ${ds.name} ${ds.n_obs}x${ds.n_vars}`)
})

await checkAsync('state: 数据集与建议', async () => {
  const state = await impl.opState(runtime)
  assert.equal(state.ok, true)
  assert.equal(state.active, datasetId)
  assert.equal(state.datasets.length, 1)
  assert.ok(Array.isArray(state.suggestions) && state.suggestions.length > 0, '应给出建议下一步')
  assert.ok(Array.isArray(state.activity) && state.activity.length > 0, '活动日志应非空')
  assert.ok(Array.isArray(state.recipes) && state.recipes.length >= 3, `预设流程不足：${JSON.stringify(state.recipes)}`)
  assert.ok(Array.isArray(state.datasets[0].history), 'datasets[].history 缺失')
  console.log(`      recipes=${state.recipes.map((r) => r.id).join(',')}`)
})

await checkAsync('audit: 只载入未分析时应给 warn/fail', async () => {
  const result = await impl.opAudit(runtime, { dataset: datasetId, actor: 'user' })
  assert.equal(result.ok, true)
  assert.ok(Number.isFinite(result.score), 'score 缺失')
  assert.ok(result.checks.length >= 8, `检查项太少：${result.checks.length}`)
  const ids = result.checks.map((c) => c.id)
  for (const id of ['qc_before_cluster', 'clusters_few', 'embedding_missing']) {
    assert.ok(ids.includes(id), `缺少检查项 ${id}`)
  }
  const problems = result.checks.filter((c) => c.level !== 'ok').length
  assert.ok(problems > 0, '空数据状态下应有问题项')
  console.log(`      score=${result.score} level=${result.level} problems=${problems}`)
})

const pipeline = ['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'embed']
for (const step of pipeline) {
  await checkAsync(`run: ${step}`, async () => {
    const result = await impl.opRun(runtime, { dataset: datasetId, step, params: {}, actor: 'user' })
    assert.equal(result.ok, true, `步骤 ${step} 失败`)
    assert.ok(result.summary && typeof result.summary.title === 'string', `${step} 缺 summary.title`)
    assert.ok(result.steps.includes(step), `${step} 未出现在 steps：${result.steps}`)
    const bullets = result.summary.bullets || []
    console.log(`      ${step}: ${bullets.slice(0, 2).join(' | ')}`)
  })
}

await checkAsync('embedding: cluster 着色坐标', async () => {
  const result = await impl.opEmbedding(runtime, { dataset: datasetId, color_by: { kind: 'cluster', key: 'clusters' } })
  assert.equal(result.ok, true)
  assert.equal(result.x.length, result.n_points, 'x 长度与 n_points 不一致')
  assert.equal(result.y.length, result.n_points, 'y 长度与 n_points 不一致')
  assert.ok(result.x.every((v) => Number.isFinite(v)), 'x 含非有限值')
  assert.ok(Array.isArray(result.categories) && result.categories.length > 1, '应有多于一个簇')
  assert.ok(result.categories.every((c) => c.color && c.label !== undefined), '分类缺 color/label')
  console.log(`      layout=${result.layout} n=${result.n_points} clusters=${result.categories.length}`)
})

await checkAsync('embedding: 按基因着色', async () => {
  const preview = await impl.opPreview(runtime, { dataset: datasetId })
  const gene = preview.genes[0]
  const result = await impl.opEmbedding(runtime, { dataset: datasetId, color_by: { kind: 'gene', key: gene } })
  assert.equal(result.ok, true)
  assert.equal(result.continuous, true, '基因着色应为连续型')
  assert.equal(result.values.length, result.n_points, 'values 长度不一致')
  console.log(`      gene=${gene} values[0]=${result.values[0]}`)
})

await checkAsync('markers: 各簇 top 基因', async () => {
  const result = await impl.opMarkers(runtime, { dataset: datasetId, top_n: 8, actor: 'user' })
  assert.equal(result.ok, true)
  assert.ok(result.groups.length > 1, '簇数为 1，聚类可能失败')
  for (const g of result.groups) {
    assert.ok(g.genes.length > 0, `簇 ${g.id} 没有标记基因`)
    const first = g.genes[0]
    assert.ok(first.gene && Number.isFinite(first.logfc), `簇 ${g.id} 标记基因字段不全`)
  }
  console.log(`      groups=${result.groups.length} top=${result.groups[0].genes[0].gene}`)
})

await checkAsync('plot: qc_violin 生成 PNG', async () => {
  const result = await impl.opPlot(runtime, { dataset: datasetId, kind: 'qc_violin', params: { metrics: ['n_genes', 'n_counts', 'pct_counts_mt'] }, actor: 'user' })
  assert.equal(result.ok, true)
  assert.ok(result.figure && result.figure.url.startsWith('/singlecell/figure?id='), 'figure.url 非法')
  assert.ok(existsSync(result.figure.file), `PNG 不存在：${result.figure.file}`)
  const magic = readFileSync(result.figure.file).subarray(0, 8)
  assert.deepEqual(Array.from(magic), [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a], 'PNG magic 不对')
  assert.ok(result.figure.width >= 200 && result.figure.height >= 200, 'PNG 尺寸过小')
  figureFile = result.figure.file
  console.log(`      figure=${path.basename(figureFile)} ${result.figure.width}x${result.figure.height}`)
})

for (const kind of ['embedding', 'dotplot', 'cluster_bar', 'pca_variance']) {
  await checkAsync(`plot: ${kind}`, async () => {
    const result = await impl.opPlot(runtime, { dataset: datasetId, kind, params: {}, actor: 'agent' })
    assert.equal(result.ok, true)
    assert.ok(existsSync(result.figure.file), `${kind} PNG 不存在`)
  })
}

await checkAsync('figure: 注册表可查', async () => {
  const fig = impl.opFigure(runtime, path.basename(figureFile))
  assert.equal(fig.file, figureFile)
  assert.throws(() => impl.opFigure(runtime, '../../etc/passwd'), /非法|不存在/, '应拒绝路径穿越')
})

await checkAsync('annotate: 写回细胞类型', async () => {
  const markers = await impl.opMarkers(runtime, { dataset: datasetId, top_n: 3 })
  const mapping = {}
  markers.groups.forEach((g, i) => { mapping[String(g.id)] = `Type-${i + 1}` })
  const result = await impl.opRun(runtime, { dataset: datasetId, step: 'annotate', params: { mapping }, actor: 'user' })
  assert.equal(result.ok, true)
  const embed = await impl.opEmbedding(runtime, { dataset: datasetId, color_by: { kind: 'obs', key: 'cell_type' } })
  assert.equal(embed.ok, true)
  assert.equal(embed.categories.length, markers.groups.length, 'cell_type 分类数与簇数不一致')
})

await checkAsync('export: h5ad + markers_csv', async () => {
  const out = path.join(workdir, 'exported.h5ad')
  const h5ad = await impl.opExport(runtime, { dataset: datasetId, what: 'h5ad', path: out, actor: 'user' })
  assert.equal(h5ad.ok, true)
  assert.ok(existsSync(out), 'h5ad 未写出')
  assert.ok(h5ad.bytes > 1000, 'h5ad 太小')
  const csv = await impl.opExport(runtime, { dataset: datasetId, what: 'markers_csv', actor: 'user' })
  assert.equal(csv.ok, true)
  assert.ok(existsSync(csv.file), 'markers csv 未写出')
  console.log(`      h5ad=${Math.round(h5ad.bytes / 1024)}KB csv=${csv.file}`)
})

await checkAsync('tool: singlecell action=status / run / markers 摘要', async () => {
  const status = await impl.runTool(runtime, { action: 'status' })
  assert.ok(status.message.includes('单细胞工作台'), 'status.message 应含工作台标题')
  assert.ok(status.active, 'status.active 缺失')
  const markers = await impl.runTool(runtime, { action: 'markers', dataset: datasetId, top_n: 3 })
  assert.ok(markers.message.includes('簇'), 'markers.message 应含簇信息')
  assert.ok(Array.isArray(markers.groups), 'markers.groups 缺失')
  console.log('      status:\n' + status.message.split('\n').map((l) => '        ' + l).join('\n'))
})

await checkAsync('prompt section: 有数据时非空、简洁', async () => {
  const text = impl.renderPromptSection(runtime)
  assert.ok(text.length > 0, '有数据集时 prompt section 不应为空')
  assert.ok(text.length < 800, `prompt section 过长：${text.length}`)
})

await checkAsync('错误路径: 未知步骤 / 未知图类型', async () => {
  await assert.rejects(() => impl.opRun(runtime, { dataset: datasetId, step: 'nope' }), /未知步骤/)
  await assert.rejects(() => impl.opPlot(runtime, { dataset: datasetId, kind: 'nope' }), /未知图类型/)
})

await checkAsync('autoannotate: 参考面板建议', async () => {
  const result = await impl.opAutoannotate(runtime, { dataset: datasetId, actor: 'user' })
  assert.equal(result.ok, true)
  assert.ok(result.suggestions.length > 1, '应给多个簇建议')
  for (const s of result.suggestions) {
    assert.ok(s.cluster !== undefined && s.label, `建议字段不全：${JSON.stringify(s)}`)
    assert.ok(Number.isFinite(s.score), `cluster ${s.cluster} score 非法`)
    assert.ok(Array.isArray(s.evidence) && s.evidence.length > 0, `cluster ${s.cluster} 缺 evidence`)
  }
  console.log('      ' + result.suggestions.slice(0, 4).map((s) => `${s.cluster}→${s.label}(${Number(s.score).toFixed(2)})`).join(' '))
})

await checkAsync('resolution_scan: 簇数与 ARI，不改当前聚类', async () => {
  const before = await impl.opMarkers(runtime, { dataset: datasetId, top_n: 2 })
  const result = await impl.opResolutionScan(runtime, { dataset: datasetId, resolutions: [0.4, 0.8, 1.2], actor: 'user' })
  assert.equal(result.ok, true)
  assert.equal(result.results.length, 3, '结果条数与输入不一致')
  assert.ok(result.current && Number.isFinite(result.current.resolution), 'current 缺失')
  for (const r of result.results) {
    assert.ok(Number.isFinite(r.n_clusters) && r.n_clusters >= 2, `簇数非法：${r.n_clusters}`)
    assert.ok(r.ari_vs_current >= -1 && r.ari_vs_current <= 1, `ARI 越界：${r.ari_vs_current}`)
  }
  const after = await impl.opMarkers(runtime, { dataset: datasetId, top_n: 2 })
  assert.equal(after.groups.length, before.groups.length, '扫描改变了当前聚类')
  console.log('      ' + result.results.map((r) => `${r.resolution}→${r.n_clusters}簇(ARI ${Number(r.ari_vs_current).toFixed(2)})`).join(' '))
})

await checkAsync('audit: 全流程后分数应更高', async () => {
  const result = await impl.opAudit(runtime, { dataset: datasetId, actor: 'user' })
  assert.equal(result.ok, true)
  const missing = result.checks.filter((c) => c.level === 'fail').map((c) => c.id)
  assert.ok(result.score >= 60, `全流程后 score 仍偏低：${result.score}（fail: ${missing.join(',')}）`)
  console.log(`      score=${result.score} level=${result.level}`)
})

await checkAsync('recipe: 未知预设报错、quick 可跑', async () => {
  await assert.rejects(() => impl.opRecipe(runtime, { dataset: datasetId, id: 'nope' }), /未知预设流程/)
  const result = await impl.opRecipe(runtime, { dataset: datasetId, id: 'quick', actor: 'user' })
  assert.equal(result.ok, true)
  assert.ok(result.completed === result.total && result.total >= 5, `quick 未跑完：${result.completed}/${result.total}`)
  console.log(`      quick ${result.completed}/${result.total} steps`)
})

await checkAsync('export: protocol / methods / report / script', async () => {
  const protocol = await impl.opExport(runtime, { dataset: datasetId, what: 'protocol', actor: 'user' })
  assert.ok(existsSync(protocol.file), 'protocol 未写出')
  const protocolJson = JSON.parse(readFileSync(protocol.file, 'utf8'))
  assert.ok(Array.isArray(protocolJson.history) && protocolJson.history.length > 0, 'protocol 缺 history')
  assert.ok(protocolJson.engine && protocolJson.engine.mode, 'protocol 缺 engine 信息')

  const methods = await impl.opExport(runtime, { dataset: datasetId, what: 'methods', actor: 'user' })
  const methodsText = readFileSync(methods.file, 'utf8')
  assert.ok(methodsText.length > 100, 'methods 太短')
  assert.ok(/cluster|Clustering/i.test(methodsText), 'methods 未提到聚类')

  const report = await impl.opExport(runtime, { dataset: datasetId, what: 'report', actor: 'user' })
  const reportText = readFileSync(report.file, 'utf8')
  assert.ok(reportText.length > 200, 'report 太短')
  assert.ok(/\.png/.test(reportText), 'report 未引用图片')

  const script = await impl.opExport(runtime, { dataset: datasetId, what: 'script', actor: 'user' })
  const scriptText = readFileSync(script.file, 'utf8')
  assert.ok(scriptText.length > 200, 'script 太短')
  assert.ok(/import /.test(scriptText), 'script 不像可运行脚本')
  console.log(`      protocol=${Math.round(protocol.bytes / 1024)}KB methods=${methodsText.length}B report=${reportText.length}B script=${scriptText.length}B`)
})

await checkAsync('state: history 溯源字段齐全', async () => {
  const state = await impl.opState(runtime)
  const ds = state.datasets[0]
  assert.ok(Array.isArray(ds.history) && ds.history.length > 0, 'history 为空')
  const last = ds.history[ds.history.length - 1]
  for (const key of ['at', 'step', 'params', 'ok', 'engine']) {
    assert.ok(key in last, `history 缺字段 ${key}`)
  }
  console.log(`      history=${ds.history.length} 条，最后一步 ${last.step}`)
})

await checkAsync('reset + dispose', async () => {
  const reset = await impl.opReset(runtime, { actor: 'user' })
  assert.equal(reset.ok, true)
  assert.equal(reset.datasets.length, 0)
  await runtime.engine.dispose()
})

const failed = results.filter((r) => !r.ok)
console.log(`\n${results.length - failed.length}/${results.length} 通过，用时 ${((Date.now() - t0) / 1000).toFixed(1)}s`)
try { rmSync(workdir, { recursive: true, force: true }) } catch { /* ignore */ }
process.exit(failed.length ? 1 : 0)
