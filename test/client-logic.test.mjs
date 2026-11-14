/**
 * dsh-singlecell-workbench 客户端纯函数单测（不依赖浏览器）。
 *
 * 运行：node test/client-logic.test.mjs
 *
 * client.js 是浏览器模块（window.__ModuleLoader__.load），但在 node 里
 * （typeof window === 'undefined'）只会把纯函数挂到 globalThis.__SC_TEST__，
 * 所以这里直接 import 它、断言纯函数行为 + 静态合规约束。
 */
import assert from 'node:assert/strict'
import { readFile } from 'node:fs/promises'
import { test } from 'node:test'

const CLIENT_URL = new URL('../client.js', import.meta.url)
await import(CLIENT_URL.href)

const SC = globalThis.__SC_TEST__
assert.ok(SC, 'node 环境应通过 globalThis.__SC_TEST__ 暴露纯函数')

const CLIENT_SOURCE = await readFile(CLIENT_URL, 'utf8')

// ---------------------------------------------------------------------------
// 0. 加载与导出形状
// ---------------------------------------------------------------------------

test('import：无 window 时不炸，且暴露 NS/API/STEPS', () => {
  assert.equal(typeof globalThis.window, 'undefined')
  assert.equal(SC.NS, 'dsh-singlecell-workbench')
  assert.equal(SC.API, '/singlecell')
  assert.equal(SC.POLL_MS, 2500)
  assert.ok(Array.isArray(SC.STEPS) && SC.STEPS.length === 10)
  assert.ok(SC.DICT_ZH.title && SC.DICT_EN.title)
  assert.ok(SC.DICT_ZH.nav && SC.DICT_ZH.engineBuiltin && SC.DICT_ZH.engineScanpy)
  assert.ok(SC.DICT_ZH.tabRecord, 'v2 记录 tab 的文案应在字典里')
  assert.ok(SC.DICT_EN.tabRecord)
})

test('字典：所有 T()/tf()/labelKey 用到的 key 都有中文词条（不出现裸 key）', () => {
  const zhStart = CLIENT_SOURCE.indexOf('var DICT_ZH = {')
  const zhEnd = CLIENT_SOURCE.indexOf('var DICT_EN = {')
  assert.ok(zhStart > 0 && zhEnd > zhStart)
  const dictBlock = CLIENT_SOURCE.slice(zhStart, zhEnd)
  const dictKeys = new Set()
  for (const match of dictBlock.matchAll(/^ {4}([A-Za-z0-9_]+):/gm)) dictKeys.add(match[1])
  assert.ok(dictKeys.size > 150, `字典词条太少：${dictKeys.size}`)

  const usedKeys = new Set()
  for (const match of CLIENT_SOURCE.matchAll(/\bT\('([A-Za-z0-9_]+)'\)/g)) usedKeys.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/tf\(T, '([A-Za-z0-9_]+)'/g)) usedKeys.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/labelKey: '([A-Za-z0-9_]+)'/g)) usedKeys.add(match[1])
  // 动态 key（step.titleKey / statusKey / engineBadge(labelKey) 等）单独补
  for (const key of ['stepQc', 'stepFilter', 'stepNormalize', 'stepHvg', 'stepPca', 'stepNeighbors', 'stepCluster', 'stepEmbed', 'stepMarkers', 'stepAnnotate']) usedKeys.add(key)
  for (const key of ['pending', 'running', 'done', 'needPrereq']) usedKeys.add(key)
  for (const key of ['engineScanpy', 'engineBuiltin', 'engineNotReady', 'engineNone']) usedKeys.add(key)

  const missing = [...usedKeys].filter((key) => !dictKeys.has(key)).sort()
  assert.deepEqual(missing, [], `字典缺 key：${missing.join(', ')}`)
})

// ---------------------------------------------------------------------------
// 1. STEPS 顺序 / 依赖 / 默认参数
// ---------------------------------------------------------------------------

test('STEPS：顺序与依赖和 CONTRACT §2.2 一致', () => {
  assert.deepEqual(SC.STEPS.map((step) => step.id), [
    'qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'embed', 'markers', 'annotate',
  ])
  const requires = {}
  for (const step of SC.STEPS) requires[step.id] = step.requires
  assert.deepEqual(requires, {
    qc: [],
    filter: ['qc'],
    normalize: ['filter'],
    hvg: ['normalize'],
    pca: ['hvg'],
    neighbors: ['pca'],
    cluster: ['neighbors'],
    embed: ['pca'],
    markers: ['cluster'],
    annotate: ['cluster'],
  })
  for (const step of SC.STEPS) {
    assert.ok(SC.STEP_TITLES_ZH[step.id], `${step.id} 应有中文标题`)
    assert.ok(step.titleKey, `${step.id} 应有 titleKey`)
  }
})

test('STEPS：默认参数与契约表一致（qc/filter/hvg/pca/cluster/embed/markers）', () => {
  assert.deepEqual(SC.defaultParams('qc'), { min_genes: 200, max_genes: null, min_counts: null, max_pct_mt: 20 })
  assert.deepEqual(SC.defaultParams('filter'), {
    min_genes: 200, max_genes: null, min_counts: null, max_pct_mt: 20, min_cells: 3,
  })
  assert.deepEqual(SC.defaultParams('normalize'), { target_sum: 10000 })
  assert.deepEqual(SC.defaultParams('hvg'), { n_top_genes: 2000, flavor: 'seurat' })
  assert.deepEqual(SC.defaultParams('pca'), { n_pcs: 30 })
  assert.deepEqual(SC.defaultParams('neighbors'), { n_neighbors: 15, metric: 'cosine' })
  assert.deepEqual(SC.defaultParams('cluster'), { resolution: 1.0, method: 'louvain' })
  assert.deepEqual(SC.defaultParams('embed'), { layout: 'tsne', perplexity: 30 })
  assert.deepEqual(SC.defaultParams('markers'), { groupby: 'clusters', top_n: 10, method: 't-test' })
  assert.deepEqual(SC.defaultParams('annotate'), {})
  assert.deepEqual(SC.defaultParams('unknown-step'), {})
})

test('collectStepParams：数字/空值/非法值', () => {
  const qc = SC.collectStepParams('qc', { min_genes: '300', max_genes: '', min_counts: null, max_pct_mt: '15' })
  assert.equal(qc.error, '')
  assert.deepEqual(qc.params, { min_genes: 300, max_genes: null, min_counts: null, max_pct_mt: 15 })

  const cluster = SC.collectStepParams('cluster', { resolution: 0.8, method: 'louvain' })
  assert.deepEqual(cluster.params, { resolution: 0.8, method: 'louvain' })

  const bad = SC.collectStepParams('pca', { n_pcs: 'abc' })
  assert.match(bad.error, /需要是数字/)
})

// ---------------------------------------------------------------------------
// 2. nextRunnableSteps / 状态
// ---------------------------------------------------------------------------

test('nextRunnableSteps：依赖满足才可运行，且按 STEPS 顺序返回', () => {
  assert.deepEqual(SC.nextRunnableSteps([]).map((step) => step.id), ['qc'])
  assert.deepEqual(SC.nextRunnableSteps(['qc']).map((step) => step.id), ['filter'])
  assert.deepEqual(
    SC.nextRunnableSteps(['qc', 'filter', 'normalize', 'hvg', 'pca']).map((step) => step.id),
    ['neighbors', 'embed'],
  )
  assert.deepEqual(
    SC.nextRunnableSteps(['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster']).map((step) => step.id),
    ['embed', 'markers', 'annotate'],
  )
  assert.deepEqual(
    SC.nextRunnableSteps(['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'embed', 'markers', 'annotate']),
    [],
  )
  assert.equal(SC.nextStep([]).id, 'qc')
  assert.equal(SC.nextStep(['qc']).id, 'filter')
  assert.equal(SC.nextStep(SC.STEPS.map((step) => step.id)), null)
})

test('stepStatusOf：running > done > pending > blocked', () => {
  assert.equal(SC.stepStatusOf('qc', [], []), 'pending')
  assert.equal(SC.stepStatusOf('filter', [], []), 'blocked')
  assert.equal(SC.stepStatusOf('filter', ['qc'], []), 'pending')
  assert.equal(SC.stepStatusOf('filter', ['qc'], ['filter']), 'running')
  assert.equal(SC.stepStatusOf('filter', ['qc', 'filter'], []), 'done')
  assert.equal(SC.stepStatusOf('filter', ['qc', 'filter'], ['filter']), 'running')
  assert.equal(SC.stepStatusKey('blocked'), 'needPrereq')
  assert.equal(SC.stepStatusKey('done'), 'done')
  assert.equal(SC.stepStatusKey('pending'), 'pending')
  assert.equal(SC.stepStatusKey('running'), 'running')
})

// ---------------------------------------------------------------------------
// 3. 数值 / 百分比 / 时间格式化
// ---------------------------------------------------------------------------

test('formatNumber：万/亿/千分位/非有限值', () => {
  assert.equal(SC.formatNumber(12345), '1.2 万')
  assert.equal(SC.formatNumber(2700), '2,700')
  assert.equal(SC.formatNumber(1234.5, 2), '1,234.5')
  assert.equal(SC.formatNumber(123456789), '1.2 亿')
  assert.equal(SC.formatNumber(0.1234, 2), '0.12')
  assert.equal(SC.formatNumber(0), '0')
  assert.equal(SC.formatNumber(-1e6), '-100 万')
  assert.equal(SC.formatNumber(null), '—')
  assert.equal(SC.formatNumber(undefined), '—')
  assert.equal(SC.formatNumber(Number.NaN), '—')
  assert.equal(SC.formatNumber(Number.POSITIVE_INFINITY), '—')
})

test('formatPct / formatPValue / formatDuration / formatClock', () => {
  assert.equal(SC.formatPct(0.9234), '92.3%')
  assert.equal(SC.formatPct(1), '100%')
  assert.equal(SC.formatPct(0), '0%')
  assert.equal(SC.formatPct(0.0005, 2), '0.05%')
  assert.equal(SC.formatPct(null), '—')

  assert.equal(SC.formatPValue(1.2e-30), '1.2e-30')
  assert.equal(SC.formatPValue(0), '0')
  assert.equal(SC.formatPValue(0.05), '0.050')
  assert.equal(SC.formatPValue(0.1234), '0.123')
  assert.equal(SC.formatPValue(-1), '—')

  assert.equal(SC.formatDuration(812), '812 ms')
  assert.equal(SC.formatDuration(1500), '1.5 s')
  assert.equal(SC.formatDuration(65000), '1.1 min')
  assert.equal(SC.formatDuration(null), '')

  const localDate = new Date(2024, 0, 2, 3, 4, 5)
  assert.equal(SC.formatClock(localDate.toISOString()), '03:04:05')
  assert.equal(SC.formatClock('not-a-date'), '')
  assert.equal(SC.formatClock(''), '')
})

test('formatError：字符串/Error/空值都收敛成短文本', () => {
  assert.equal(SC.formatError(new Error('引擎未就绪')), '引擎未就绪')
  assert.equal(SC.formatError('失败了\n换行'), '失败了 换行')
  assert.equal(SC.formatError(null), '')
  assert.equal(SC.formatError(undefined), '')
  assert.ok(SC.formatError(new Error('x'.repeat(2000))).length <= 600)
})

// ---------------------------------------------------------------------------
// 4. 坐标变换（缩放/平移/适配/命中）
// ---------------------------------------------------------------------------

const CENTER = { x: 320, y: 210 }
const VIEW = { scale: 2.5, offsetX: 40, offsetY: -30 }

test('坐标变换：worldToScreen ↔ screenToWorld 可往返', () => {
  const points = [{ x: 0, y: 0 }, { x: 3.25, y: -2.5 }, { x: -100, y: 80 }]
  for (const point of points) {
    const screen = SC.worldToScreen(point, VIEW, CENTER)
    const back = SC.screenToWorld(screen, VIEW, CENTER)
    assert.ok(Math.abs(back.x - point.x) < 1e-9, `x 往返失败：${point.x} → ${back.x}`)
    assert.ok(Math.abs(back.y - point.y) < 1e-9, `y 往返失败：${point.y} → ${back.y}`)
  }
})

test('坐标变换：zoomViewAt 保持锚点下方的数据坐标不动', () => {
  const anchor = { x: 120, y: 90 }
  const worldBefore = SC.screenToWorld(anchor, VIEW, CENTER)
  const zoomed = SC.zoomViewAt(VIEW, anchor, 1.75, CENTER)
  assert.ok(Math.abs(zoomed.scale - VIEW.scale * 1.75) < 1e-9)
  const screenAfter = SC.worldToScreen(worldBefore, zoomed, CENTER)
  assert.ok(Math.abs(screenAfter.x - anchor.x) < 1e-9)
  assert.ok(Math.abs(screenAfter.y - anchor.y) < 1e-9)
  // 缩放上下限
  const tiny = SC.zoomViewAt(VIEW, anchor, 1e-9, CENTER, 0.02, 500)
  const huge = SC.zoomViewAt(VIEW, anchor, 1e9, CENTER, 0.02, 500)
  assert.equal(tiny.scale, 0.02)
  assert.equal(huge.scale, 500)
})

test('坐标变换：panView 是纯平移', () => {
  const panned = SC.panView(VIEW, 12, -8)
  assert.deepEqual(panned, { scale: 2.5, offsetX: 52, offsetY: -38 })
  const projected = SC.worldToScreen({ x: 1, y: 1 }, VIEW, CENTER)
  const projectedPanned = SC.worldToScreen({ x: 1, y: 1 }, panned, CENTER)
  assert.ok(Math.abs(projectedPanned.x - projected.x - 12) < 1e-9)
  assert.ok(Math.abs(projectedPanned.y - projected.y + 8) < 1e-9)
})

test('fitView：数据整体落在画布内且居中', () => {
  const bounds = { minX: -5, maxX: 10, minY: -2, maxY: 6 }
  const width = 800
  const height = 400
  const center = SC.viewCenter(width, height)
  const fitted = SC.fitView(bounds, width, height, 20)
  const corners = [
    { x: bounds.minX, y: bounds.minY },
    { x: bounds.maxX, y: bounds.maxY },
    { x: bounds.minX, y: bounds.maxY },
    { x: bounds.maxX, y: bounds.minY },
  ]
  for (const corner of corners) {
    const screen = SC.worldToScreen(corner, fitted, center)
    assert.ok(screen.x >= -1e-6 && screen.x <= width + 1e-6, `x 越界：${screen.x}`)
    assert.ok(screen.y >= -1e-6 && screen.y <= height + 1e-6, `y 越界：${screen.y}`)
  }
  const dataCenter = SC.worldToScreen({ x: 2.5, y: 2 }, fitted, center)
  assert.ok(Math.abs(dataCenter.x - width / 2) < 1e-6)
  assert.ok(Math.abs(dataCenter.y - height / 2) < 1e-6)
  // 退化数据（所有点相同）不能 NaN
  const degenerate = SC.fitView({ minX: 1, maxX: 1, minY: 2, maxY: 2 }, 300, 200, 10)
  assert.ok(Number.isFinite(degenerate.scale) && degenerate.scale > 0)
  assert.deepEqual(SC.fitView(null, 300, 200, 10), SC.createView())
})

test('computeBounds / findNearestPoint', () => {
  assert.deepEqual(SC.computeBounds([1, 2, 3], [4, 5, 6]), { minX: 1, maxX: 3, minY: 4, maxY: 6 })
  assert.equal(SC.computeBounds([], []), null)
  assert.equal(SC.computeBounds(['x'], ['y']), null)

  const xs = [0, 10, 20]
  const ys = [0, 10, 20]
  const view = { scale: 10, offsetX: 0, offsetY: 0 }
  const center = { x: 100, y: 100 }
  const second = SC.worldToScreen({ x: 10, y: 10 }, view, center)
  assert.equal(SC.findNearestPoint(xs, ys, view, center, { x: second.x + 2, y: second.y - 1 }, 6), 1)
  assert.equal(SC.findNearestPoint(xs, ys, view, center, { x: second.x + 400, y: second.y }, 6), -1)
  assert.equal(SC.findNearestPoint([], [], view, center, { x: 0, y: 0 }, 6), -1)
})

// ---------------------------------------------------------------------------
// 5. tooltip 文本
// ---------------------------------------------------------------------------

test('tooltipText：分类给分组名、连续给指标+数值（细胞序号从 1 开始）', () => {
  assert.equal(SC.tooltipText({ index: 12, label: '0' }), '细胞 #13 · 0')
  assert.equal(SC.tooltipText({ index: 0, label: 'T cell' }), '细胞 #1 · T cell')
  assert.equal(SC.tooltipText({ index: 3, label: '' }), '细胞 #4')
  assert.equal(SC.tooltipText({ index: 0, continuous: true, title: 'n_genes', value: 1234.5 }), '细胞 #1 · n_genes 1,234.5')
  assert.equal(SC.tooltipText({ index: 1, continuous: true, value: 0.5 }), '细胞 #2 · 0.5')
  assert.equal(SC.tooltipText(null), '细胞 #?')
})

// ---------------------------------------------------------------------------
// 6. 建议 → 中文指令 / 记录展示
// ---------------------------------------------------------------------------

test('buildAgentInstruction：带 step 的指令带参数，无 step 的退回标题/动作', () => {
  assert.equal(
    SC.buildAgentInstruction({ step: 'cluster', params: { resolution: 1 } }),
    '用单细胞工作台把当前数据集跑到聚类，resolution=1。',
  )
  assert.equal(
    SC.buildAgentInstruction({ step: 'markers', params: { top_n: 10, method: 't-test' } }),
    '用单细胞工作台把当前数据集跑到标记基因，top_n=10，method=t-test。',
  )
  assert.equal(
    SC.buildAgentInstruction({ step: 'qc', params: {} }),
    '用单细胞工作台把当前数据集跑到质控。',
  )
  assert.equal(SC.buildAgentInstruction({ id: 'demo', title: '载入示例数据', action: 'demo' }), '用单细胞工作台载入内置示例数据。')
  assert.equal(SC.buildAgentInstruction({ id: 'open', action: 'open' }), '用单细胞工作台打开一个 h5ad 数据集。')
  assert.equal(SC.buildAgentInstruction({ title: '导出结果' }), '用单细胞工作台执行「导出结果」。')
  assert.equal(SC.buildAgentInstruction(null), '用单细胞工作台继续当前数据集的分析。')
})

test('formatParamList / formatParamValue', () => {
  assert.equal(SC.formatParamList({ resolution: 1, method: 'louvain' }), 'resolution=1，method=louvain')
  assert.equal(SC.formatParamList({ mapping: { 0: 'T cell' } }), 'mapping={"0":"T cell"}')
  assert.equal(SC.formatParamList({ x: null, y: [1, 2] }), 'x=null，y=1/2')
  assert.equal(SC.formatParamList(null), '')
  assert.equal(SC.formatParamValue(true), 'true')
})

test('formatHistoryEntry：时间/步骤/参数/摘要/耗时/引擎', () => {
  const view = SC.formatHistoryEntry({
    at: new Date(2024, 0, 2, 3, 4, 5).toISOString(),
    step: 'cluster',
    params: { resolution: 1 },
    ok: true,
    ms: 812,
    summary: 'resolution=1 → 8 簇',
    n_obs: 2700,
    n_vars: 1200,
    engine: { mode: 'builtin', python: '3.14.2' },
  })
  assert.equal(view.time, '03:04:05')
  assert.equal(view.step, 'cluster')
  assert.equal(view.title, '聚类')
  assert.equal(view.paramsText, 'resolution=1')
  assert.equal(view.summary, 'resolution=1 → 8 簇')
  assert.equal(view.duration, '812 ms')
  assert.equal(view.engineMode, 'builtin')
  assert.equal(view.ok, true)
  const failed = SC.formatHistoryEntry({ step: 'pca', ok: false })
  assert.equal(failed.ok, false)
  assert.equal(failed.title, 'PCA')
})

test('buildAnnotationMapping / parseResolutions / recipeProgressText', () => {
  assert.deepEqual(
    SC.buildAnnotationMapping([
      { cluster: '0', label: 'T cell' },
      { cluster: '1', label: ' B cell ' },
      { cluster: '', label: 'skip' },
      { cluster: '2', label: '' },
      { cluster: '0', label: 'duplicate-ignored' },
    ]),
    { 0: 'T cell', 1: 'B cell' },
  )
  assert.deepEqual(SC.buildAnnotationMapping(null), {})

  assert.deepEqual(SC.parseResolutions('0.4,0.8,1.2'), { list: [0.4, 0.8, 1.2], error: '' })
  assert.deepEqual(SC.parseResolutions('0.4，0.8 1.2'), { list: [0.4, 0.8, 1.2], error: '' })
  assert.deepEqual(SC.parseResolutions('0.4,0.4,9'), { list: [0.4], error: '' })
  assert.equal(SC.parseResolutions('abc').error.length > 0, true)
  assert.equal(SC.parseResolutions([0.5, 1, 2, 3, 4, 5, 0.6, 1.5, 2.5]).list.length, 8)

  assert.equal(SC.recipeProgressText({ index: 2, total: 6, step: 'cluster' }), '第 2/6 步：聚类')
  assert.equal(SC.recipeProgressText(null), '')
})

// ---------------------------------------------------------------------------
// 7. 归一化：缺字段/旧返回不能崩
// ---------------------------------------------------------------------------

test('normalizeState：空值安全 + v2 字段（history/recipes/recipeProgress）', () => {
  const empty = SC.normalizeState(null)
  assert.equal(empty.active, null)
  assert.deepEqual(empty.datasets, [])
  assert.equal(empty.engine.mode, 'unknown')
  assert.deepEqual(empty.suggestions, [])
  assert.deepEqual(empty.recipes, [])
  assert.equal(empty.recipeProgress, null)

  const payload = {
    ok: true,
    active: 'ds1',
    engine: { ready: true, mode: 'builtin', python: '3.14.2', layouts: ['tsne'] },
    datasets: [{
      id: 'ds1', name: 'demo-pbmc', n_obs: 2700, n_vars: 1200, steps: ['qc', 'cluster'],
      current: 'cluster', warnings: ['w'],
      history: [{ at: new Date(2024, 0, 2, 3, 4, 5).toISOString(), step: 'qc', params: {}, ok: true, ms: 10, summary: 'qc ok' }],
    }],
    suggestions: [{ id: 'run:embed', title: '嵌入', step: 'embed', params: { layout: 'tsne' } }],
    recipes: [{ id: 'standard', title: '标准流程', detail: 'd', steps: [{ step: 'qc', params: {} }, { step: 'nope' }] }],
    recipeProgress: { recipe: 'standard', title: '标准流程', index: 1, total: 6, step: 'qc' },
    activity: [{ at: 'x', actor: 'agent', action: 'run:qc', detail: 'd', ok: true, ms: 10 }],
  }
  const state = SC.normalizeState(payload)
  assert.equal(state.active, 'ds1')
  assert.equal(state.datasets.length, 1)
  assert.equal(state.datasets[0].nObs, 2700)
  assert.deepEqual(state.datasets[0].steps, ['qc', 'cluster'])
  assert.equal(state.datasets[0].history.length, 1)
  assert.equal(state.engine.mode, 'builtin')
  assert.equal(state.engine.ready, true)
  assert.equal(state.suggestions[0].step, 'embed')
  assert.equal(state.recipes[0].steps.length, 1, '未知 step 应被丢弃')
  assert.equal(state.recipeProgress.index, 1)
  assert.equal(state.activity[0].actor, 'agent')
  assert.equal(SC.findDataset(state, 'ds1').id, 'ds1')
  assert.equal(SC.findDataset(state, 'nope'), null)
})

test('normalizeEmbedding：长度不一致时截断、非有限坐标归零、长度自洽', () => {
  const data = SC.normalizeEmbedding({
    x: [1, 2, 3, 'oops'],
    y: [1, 2],
    group_index: [0, 1],
    categories: [{ id: '0', label: '0', color: '#111111', n: 1 }, { id: '1', label: '1', n: 1 }],
    continuous: false,
    layout: 'tsne',
  })
  assert.equal(data.nPoints, 2)
  assert.equal(data.x.length, 2)
  assert.equal(data.y.length, 2)
  assert.equal(data.groupIndex.length, 2)
  assert.equal(data.categories[1].color, SC.PALETTE[1], '缺 color 时用调色板兜底')

  const continuous = SC.normalizeEmbedding({
    x: [1, 2, 3], y: [1, 2, 3], values: [0.5, 0.7], continuous: true, legend_title: 'CD3D',
  })
  assert.equal(continuous.nPoints, 2)
  assert.equal(continuous.values.length, 2)
  assert.equal(continuous.legendTitle, 'CD3D')
  assert.equal(SC.normalizeEmbedding(null).nPoints, 0)
})

test('normalizeMarkersGroups / normalizeFigure', () => {
  const groups = SC.normalizeMarkersGroups({ groups: [{ id: 0, label: '0', n_cells: 312, genes: [{ gene: 'CD3D', logfc: 2.3, pval_adj: 1e-30, pct_in: 0.92 }] }] })
  assert.equal(groups.length, 1)
  assert.equal(groups[0].id, '0')
  assert.equal(groups[0].nCells, 312)
  assert.equal(groups[0].genes[0].logfc, 2.3)
  assert.equal(groups[0].genes[0].pvalAdj, 1e-30)
  assert.equal(groups[0].genes[0].pctIn, 0.92)
  assert.deepEqual(SC.normalizeMarkersGroups(null), [])
  assert.deepEqual(SC.normalizeMarkersGroups({ markers: [{ id: '1', genes: [] }] }).length, 1)

  const figure = SC.normalizeFigure({ ok: true, figure: { id: 'f1', url: '/singlecell/figure?id=f1.png', width: 720, height: 480, kind: 'qc_violin' } })
  assert.equal(figure.url, '/singlecell/figure?id=f1.png')
  assert.equal(figure.width, 720)
  assert.equal(figure.kind, 'qc_violin')
  assert.equal(SC.normalizeFigure({ ok: true }), null)
})

test('v2 归一化：audit / autoannotate / resolution_scan / export', () => {
  const audit = SC.normalizeAudit({
    ok: true, score: 78.6, level: 'warn',
    checks: [{ id: 'cells_low', level: 'fail', title: '细胞过少', detail: 'd', hint: 'h' }, { id: 'unknown', level: 'weird', title: 't' }],
  })
  assert.equal(audit.score, 79)
  assert.equal(audit.level, 'warn')
  assert.equal(audit.checks[0].level, 'fail')
  assert.equal(audit.checks[1].level, 'ok', '未知等级回落 ok')
  assert.equal(audit.checks[1].detail, '')

  const autoann = SC.normalizeAutoannotate({
    panel: 'builtin-marker-panel-v1',
    suggestions: [{ cluster: '0', label: 'T cell', score: 0.82, margin: 0.19, evidence: ['CD3D'], alternative: 'NK cell' }],
  })
  assert.equal(autoann.panel, 'builtin-marker-panel-v1')
  assert.equal(autoann.suggestions[0].cluster, '0')
  assert.deepEqual(autoann.suggestions[0].evidence, ['CD3D'])
  assert.deepEqual(SC.normalizeAutoannotate(null).suggestions, [])

  const scan = SC.normalizeResolutionScan({
    current: { resolution: 1.0, n_clusters: 8 },
    results: [{ resolution: 0.4, n_clusters: 5, ari_vs_current: 0.61 }, { resolution: 0.8, n_clusters: 7 }],
  })
  assert.deepEqual(scan.current, { resolution: 1.0, nClusters: 8 })
  assert.equal(scan.results[0].ari, 0.61)
  assert.equal(scan.results[1].ari, null, '缺 ARI 时置 null（UI 显示 —）')
  assert.equal(SC.normalizeResolutionScan(null).current, null)

  const exported = SC.normalizeExportResult({ ok: true, what: 'report', file: '/tmp/x/report.md', bytes: 1234 })
  assert.equal(exported.file, '/tmp/x/report.md')
  assert.equal(exported.bytes, 1234)
  assert.equal(SC.normalizeExportResult(null).file, '')
})

test('normalizePreview / normalizeActivity / engineBadge', () => {
  const preview = SC.normalizePreview({ genes: ['CD3D'], obs_columns: [{ key: 'n_genes', kind: 'numeric' }], n_obs: 10, n_vars: 20, steps: ['qc'] })
  assert.deepEqual(preview.genes, ['CD3D'])
  assert.equal(preview.obsColumns[0].kind, 'numeric')
  assert.deepEqual(SC.normalizePreview(null).genes, [])

  const activity = SC.normalizeActivity([{ at: 'x', actor: 'agent', action: 'run:cluster', detail: 'd', ok: false, ms: 12 }])
  assert.equal(activity[0].actor, 'agent')
  assert.equal(activity[0].ok, false)
  assert.equal(activity[0].ms, 12)
  assert.equal(SC.normalizeActivity(null).length, 0)

  assert.equal(SC.engineBadge({ mode: 'scanpy', python: '3.12.0' }).labelKey, 'engineScanpy')
  assert.equal(SC.engineBadge({ mode: 'builtin', python: '3.14.2' }).tone, 'ok')
  assert.equal(SC.engineBadge({ mode: 'none', error: '没装 numpy' }).tone, 'bad')
  assert.equal(SC.engineBadge({ mode: 'unknown', ready: false }).labelKey, 'engineNotReady')

  const chosen = SC.chooseEngine(SC.normalizeEngine({ mode: 'unknown' }), SC.normalizeEngine({ ready: true, mode: 'builtin' }))
  assert.equal(chosen.mode, 'builtin')
})

test('v2 常量：AUDIT_LEVEL_STYLE / EXPORT_KINDS', () => {
  assert.equal(SC.AUDIT_LEVEL_STYLE.ok.tone, 'ok')
  assert.equal(SC.AUDIT_LEVEL_STYLE.warn.tone, 'warn')
  assert.equal(SC.AUDIT_LEVEL_STYLE.fail.tone, 'fail')
  assert.equal(SC.AUDIT_LEVEL_STYLE.warn.icon, '!')
  assert.deepEqual(SC.EXPORT_KINDS.map((kind) => kind.id), ['protocol', 'methods', 'report', 'script', 'h5ad', 'markers_csv'])
})

test('pickContinuousColor：返回 rgb 且随值单调变化', () => {
  const low = SC.pickContinuousColor(0, 0, 1)
  const mid = SC.pickContinuousColor(0.5, 0, 1)
  const high = SC.pickContinuousColor(1, 0, 1)
  for (const color of [low, mid, high]) assert.match(color, /^rgb\(\d+,\d+,\d+\)$/)
  assert.equal(low, 'rgb(59,76,192)')
  assert.equal(high, 'rgb(192,57,43)')
  assert.notEqual(low, mid)
  assert.notEqual(mid, high)
  assert.equal(SC.pickContinuousColor(Number.NaN, 0, 1), SC.CONTINUOUS_STOPS[0])
})

// ---------------------------------------------------------------------------
// 8. API 客户端：假 fetch 下的 method/url/body 与错误解析
// ---------------------------------------------------------------------------

function fakeResponse(body, { status = 200, ok = status < 400 } = {}) {
  const text = typeof body === 'string' ? body : JSON.stringify(body)
  return { ok, status, text: async () => text }
}

function createFakeFetch(responses) {
  const calls = []
  const queue = Array.isArray(responses) ? responses.slice() : []
  const fake = async (url, options) => {
    calls.push({ url, options })
    const next = queue.length ? queue.shift() : { body: { ok: true } }
    if (next instanceof Error) throw next
    return fakeResponse(next.body, next)
  }
  return { fake, calls }
}

test('api：所有请求都走同源相对路径 /singlecell/*', async () => {
  const { fake, calls } = createFakeFetch([])
  const api = SC.createApi(fake)
  assert.equal(api.base, '/singlecell')

  await api.getHealth()
  await api.getState()
  await api.getActivity()
  await api.openPath('/tmp/a.h5ad')
  await api.openDemo()
  await api.selectDataset('ds1')
  await api.closeDataset('ds1')
  await api.preview('ds1')
  await api.runStep('ds1', 'cluster', { resolution: 1 })
  await api.getEmbedding('ds1', { kind: 'cluster', key: 'clusters' })
  await api.plot('ds1', 'qc_violin', { metrics: ['n_genes'] })
  await api.markers('ds1', { top_n: 10, method: 't-test' })
  await api.getAudit('ds1')
  await api.getAutoannotate('ds1')
  await api.runRecipe('ds1', 'standard')
  await api.resolutionScan('ds1', [0.4, 0.8])
  await api.exportData('ds1', 'report')
  await api.reset()

  assert.equal(calls.length, 18)
  for (const call of calls) {
    assert.ok(call.url.startsWith('/singlecell/'), `非相对路径：${call.url}`)
    assert.ok(!/^https?:/.test(call.url), `不应出现绝对 URL：${call.url}`)
  }
  assert.deepEqual(calls.map((call) => call.url), [
    '/singlecell/health',
    '/singlecell/state',
    '/singlecell/activity',
    '/singlecell/dataset/open',
    '/singlecell/dataset/demo',
    '/singlecell/dataset/select',
    '/singlecell/dataset/close',
    '/singlecell/preview',
    '/singlecell/run',
    '/singlecell/embedding',
    '/singlecell/plot',
    '/singlecell/markers',
    '/singlecell/audit',
    '/singlecell/autoannotate',
    '/singlecell/recipe',
    '/singlecell/resolution-scan',
    '/singlecell/export',
    '/singlecell/reset',
  ])
})

test('api GET：method=GET，无 body', async () => {
  const { fake, calls } = createFakeFetch([{ body: { ok: true, active: 'ds1', datasets: [] } }])
  const api = SC.createApi(fake)
  const payload = await api.getState()
  assert.equal(payload.active, 'ds1')
  assert.equal(calls[0].url, '/singlecell/state')
  assert.equal(calls[0].options.method, 'GET')
  assert.equal(calls[0].options.body, undefined)
})

test('api POST /run：JSON body 含 dataset/step/params/actor', async () => {
  const { fake, calls } = createFakeFetch([{ body: { ok: true, step: 'cluster', steps: ['qc'] } }])
  const api = SC.createApi(fake)
  await api.runStep('ds1', 'cluster', { resolution: 0.8 })
  const options = calls[0].options
  assert.equal(calls[0].url, '/singlecell/run')
  assert.equal(options.method, 'POST')
  assert.equal(options.headers['content-type'], 'application/json')
  assert.deepEqual(JSON.parse(options.body), {
    dataset: 'ds1',
    step: 'cluster',
    params: { resolution: 0.8 },
    actor: 'user',
  })
})

test('api v2 路由：audit / autoannotate / recipe / resolution-scan / export body 正确', async () => {
  const { fake, calls } = createFakeFetch([])
  const api = SC.createApi(fake)
  await api.getAudit('ds1')
  await api.getAutoannotate('ds1')
  await api.runRecipe('ds1', 'standard')
  await api.resolutionScan('ds1', [0.4, 1.2])
  await api.exportData('ds1', 'methods')
  await api.getEmbedding('ds1', { kind: 'obs', key: 'n_genes' })

  assert.deepEqual(JSON.parse(calls[0].options.body), { dataset: 'ds1', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[1].options.body), { dataset: 'ds1', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[2].options.body), { dataset: 'ds1', id: 'standard', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[3].options.body), { dataset: 'ds1', resolutions: [0.4, 1.2], actor: 'user' })
  assert.deepEqual(JSON.parse(calls[4].options.body), { dataset: 'ds1', what: 'methods', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[5].options.body), { dataset: 'ds1', color_by: { kind: 'obs', key: 'n_genes' } })
})

test('api upload：raw File body + x-dsh-filename header，不设 content-type', async () => {
  const file = { name: 'sample.h5ad', size: 1024 }
  const { fake, calls } = createFakeFetch([{ body: { ok: true, dataset: { id: 'ds9' } } }])
  const api = SC.createApi(fake)
  await api.upload(file)
  assert.equal(calls[0].url, '/singlecell/dataset/upload')
  assert.equal(calls[0].options.method, 'POST')
  assert.equal(calls[0].options.body, file)
  assert.equal(calls[0].options.headers['x-dsh-filename'], 'sample.h5ad')
  assert.equal(calls[0].options.headers['content-type'], undefined)
})

test('api 错误：{ok:false,error} / HTTP 错误 / 无 error 字段 / 网络异常都转成中文 Error', async () => {
  const prerequisite = createFakeFetch([{ body: { ok: false, code: 'prerequisite', error: '缺少前置步骤 qc' } }])
  await assert.rejects(
    () => SC.createApi(prerequisite.fake).getEmbedding('ds1', { kind: 'cluster', key: 'clusters' }),
    /缺少前置步骤 qc/,
  )

  const httpError = createFakeFetch([{ body: { ok: false, error: '引擎未就绪' }, status: 500 }])
  await assert.rejects(() => SC.createApi(httpError.fake).getState(), /引擎未就绪/)

  const noDetail = createFakeFetch([{ body: { ok: false, code: 'internal' } }])
  await assert.rejects(() => SC.createApi(noDetail.fake).getState(), /请求失败/)

  const network = createFakeFetch([new TypeError('fetch failed')])
  await assert.rejects(() => SC.createApi(network.fake).getState(), /网络请求失败/)

  // 返回非 JSON 文本 + 500 时至少给 HTTP 状态
  const plainText = createFakeFetch([{ body: 'boom', status: 502 }])
  await assert.rejects(() => SC.createApi(plainText.fake).getState(), /HTTP 502/)
})

test('api 超时：可读中文错误（注入 10ms 超时 + 响应 abort）', async () => {
  let aborted = false
  const fake = (url, options) => new Promise((resolve, reject) => {
    const signal = options && options.signal
    if (signal && typeof signal.addEventListener === 'function') {
      signal.addEventListener('abort', () => {
        aborted = true
        const error = new Error('Aborted')
        error.name = 'AbortError'
        reject(error)
      })
    }
  })
  const api = SC.createApi(fake)
  await assert.rejects(() => api.getState({ timeoutMs: 10 }), /请求超时/)
  assert.equal(aborted, true)
})

// ---------------------------------------------------------------------------
// 8.5 v3：全格式导入 + 高级分析（纯函数）
// ---------------------------------------------------------------------------

test('ANALYSIS_GROUPS：10 个分组且顺序与 registry/GROUPS 一致', () => {
  assert.deepEqual(SC.ANALYSIS_GROUPS, [
    '基本分析', '组间对比', '数据可视化', '通路分析', '基因集合分析',
    '发育轨迹分析', '肿瘤异质性分析', '细胞通讯分析', '转录调控分析', '虚拟敲除',
  ])
  assert.equal(SC.ANALYSIS_GROUPS.length, 10)
  for (const group of SC.ANALYSIS_GROUPS) {
    assert.ok(SC.ANALYSIS_GROUP_LABEL_KEYS[group], `${group} 应有 groupXxx 字典 key`)
  }
  assert.equal(SC.ANALYSIS_GROUP_OTHER, '其他')
})

test('normalizeAnalyses：native/external、available/reason、params/needs 归一化', () => {
  const items = SC.normalizeAnalyses({
    analyses: [
      {
        id: 'cell_cycle', title: '细胞周期评分', group: '基本分析', method: 'Tirosh 打分',
        native: true, available: true, needs: ['clusters'],
        params: [{ key: 'n', label: '数量', type: 'number', default: 6, min: 2, max: 20, hint: 'h' }],
      },
      { id: 'de_group', title: '组间差异', group: '组间对比', native: true, available: false, reason: '需要 group 列' },
      { id: 'rna_velocity', title: 'RNA velocity', group: '发育轨迹分析', native: false, reason: '需要剪接层', alternative: '用 trajectory_dpt' },
      { title: '没有 id 的条目' },
    ],
  })
  assert.equal(items.length, 3, '缺 id 的条目应被丢弃')
  assert.equal(items[0].id, 'cell_cycle')
  assert.equal(items[0].native, true)
  assert.equal(items[0].external, false)
  assert.equal(items[0].available, true)
  assert.deepEqual(items[0].needs, ['clusters'])
  assert.equal(items[0].params[0].key, 'n')
  assert.equal(items[0].params[0].type, 'number')
  assert.equal(items[0].params[0].min, 2)
  assert.equal(items[1].available, false)
  assert.equal(items[1].reason, '需要 group 列')
  assert.equal(items[2].native, false)
  assert.equal(items[2].external, true)
  assert.equal(items[2].available, false)
  assert.equal(items[2].alternative, '用 trajectory_dpt')
  assert.deepEqual(SC.normalizeAnalyses(null), [])
  assert.deepEqual(SC.normalizeAnalyses([{ id: 'x', native: true }])[0].params, [])
})

test('groupAnalyses：按固定分组顺序，未知分组排最后', () => {
  const items = [
    { id: 'a', group: '虚拟敲除' },
    { id: 'b', group: '基本分析' },
    { id: 'c', group: '自定义组' },
    { id: 'd', group: '基本分析' },
    { id: 'e' },
  ]
  const groups = SC.groupAnalyses(items)
  assert.deepEqual(groups.map((entry) => entry.group), ['基本分析', '虚拟敲除', '自定义组', SC.ANALYSIS_GROUP_OTHER])
  assert.deepEqual(groups[0].items.map((item) => item.id), ['b', 'd'])
  assert.equal(groups[0].labelKey, 'groupBasic')
  assert.equal(groups[3].labelKey, '')
  assert.deepEqual(SC.groupAnalyses(null), [])
})

test('paramDefaults / paramOptionValues：number/boolean/select/column/file/string 默认值', () => {
  const params = [
    { key: 'n', type: 'number', default: 5 },
    { key: 'n2', type: 'number' },
    { key: 'b', type: 'boolean', default: true },
    { key: 's', type: 'select', options: ['x', 'y'], default: 'y' },
    { key: 's2', type: 'select', options: [{ value: 'a', label: 'A' }] },
    { key: 'c', type: 'column', default: 'clusters' },
    { key: 'f', type: 'file', default: '/tmp/x.gmt' },
    { key: 't', type: 'string', default: 'hi' },
  ]
  assert.deepEqual(SC.paramDefaults(params), {
    n: 5, n2: '', b: true, s: 'y', s2: 'a', c: 'clusters', f: '/tmp/x.gmt', t: 'hi',
  })
  assert.deepEqual(SC.paramOptionValues({ type: 'select', options: ['x', { value: 'a', label: 'A' }] }), [
    { value: 'x', label: 'x' },
    { value: 'a', label: 'A' },
  ])
  assert.deepEqual(SC.paramDefaults(null), {})
})

test('collectAnalysisParams：空值省略、boolean 归一、非法数字给可读错误', () => {
  const params = [
    { key: 'n', type: 'number' },
    { key: 'n2', type: 'number' },
    { key: 'b', type: 'boolean' },
    { key: 's', type: 'select' },
    { key: 'c', type: 'column' },
    { key: 'f', type: 'file' },
  ]
  const built = SC.collectAnalysisParams(params, { n: '7', n2: '', b: false, s: 'y', c: 'clusters', f: '' })
  assert.equal(built.error, '')
  assert.deepEqual(built.params, { n: 7, b: false, s: 'y', c: 'clusters' })
  const bad = SC.collectAnalysisParams(params, { n: 'abc' })
  assert.match(bad.error, /需要是数字/)
})

test('needsHint：clusters / groups:<col> / gene_order / 空 needs', () => {
  const noSteps = { steps: [], obsColumns: [] }
  assert.equal(SC.needsHint({ needs: [] }, noSteps), '')
  assert.match(SC.needsHint({ needs: ['clusters'] }, noSteps), /需要先运行聚类步骤/)
  assert.equal(SC.needsHint({ needs: ['clusters'] }, { steps: ['cluster'], obsColumns: [] }), '')
  assert.match(SC.needsHint({ needs: ['groups:group'] }, { steps: ['cluster'], obsColumns: [{ key: 'clusters' }] }), /group/)
  assert.equal(SC.needsHint({ needs: ['groups:group'] }, { steps: ['cluster'], obsColumns: [{ key: 'group', kind: 'categorical' }] }), '')
  assert.match(SC.needsHint({ needs: ['gene_order'] }, noSteps), /基因位置/)
  assert.equal(SC.needsHint(null, null), '')
})

test('formatCandidates：archive member / RData object / 字符串候选 / 嵌套 data.candidates', () => {
  const archive = SC.formatCandidates({
    ok: false,
    code: 'ambiguous',
    error: '多个候选',
    candidates: [
      { member: 'counts/matrix.mtx', label: 'matrix.mtx', size: 123456, format: '10x_mtx' },
      { member: 'GSE1_series_matrix.txt.gz', size: 98765 },
    ],
  })
  assert.equal(archive.count, 2)
  assert.equal(archive.items[0].label, 'matrix.mtx')
  assert.deepEqual(archive.items[0].patch, { member: 'counts/matrix.mtx' })
  assert.match(archive.items[0].detail, /10x_mtx/)
  assert.match(archive.items[0].detail, /KB/)
  assert.deepEqual(archive.items[1].patch, { member: 'GSE1_series_matrix.txt.gz' })

  const rdata = SC.formatCandidates({
    candidates: [{ object: 'seu', class: 'Seurat', label: 'seu', n_obs: 2400, n_vars: 1200, size: 2048 }],
  })
  assert.deepEqual(rdata.items[0].patch, { object: 'seu' })
  assert.match(rdata.items[0].detail, /Seurat · 2,400×1,200 · 2 KB/)

  const nested = SC.formatCandidates({ data: { candidates: [{ member: 'a/b.csv' }] } })
  assert.deepEqual(nested.items[0].patch, { member: 'a/b.csv' })

  const strings = SC.formatCandidates(['a/b.tsv', 'someObject'])
  assert.equal(strings.items.length, 2)
  assert.deepEqual(strings.items[0].patch, { member: 'a/b.tsv' })
  assert.equal(SC.formatCandidates(null).count, 0)
})

test('formatTableRows / formatTableCell：列、单元格展示值、行数与截断', () => {
  const table = {
    id: 't1', title: 'DE', note: 'n',
    columns: [
      { key: 'gene', label: 'Gene', type: 'string' },
      { key: 'logfc', label: 'logFC', type: 'number' },
    ],
    rows: [{ gene: 'CD3D', logfc: 2.3456 }, { gene: 'CD3E', logfc: 1.2 }],
  }
  const view = SC.formatTableRows(table)
  assert.equal(view.total, 2)
  assert.equal(view.shown, 2)
  assert.equal(view.truncated, false)
  assert.equal(view.columns[1].numeric, true)
  assert.equal(view.rows[0].cells[0].text, 'CD3D')
  assert.equal(view.rows[0].cells[1].text, '2.346')
  const limited = SC.formatTableRows(table, 1)
  assert.equal(limited.shown, 1)
  assert.equal(limited.truncated, true)
  const broken = SC.formatTableRows(null)
  assert.deepEqual([broken.columns.length, broken.rows.length, broken.total, broken.truncated], [0, 0, 0, false])

  assert.equal(SC.formatTableCell(null, true), '—')
  assert.equal(SC.formatTableCell(2.3456, true), '2.346')
  assert.equal(SC.formatTableCell(2.5, false), '2.5')
  assert.equal(SC.formatTableCell(true, true), 'true')
  assert.equal(SC.formatTableCell(['a', 'b'], true), 'a, b')
})

test('formatBytes：B/KB/MB 与非法值', () => {
  assert.equal(SC.formatBytes(0), '0 B')
  assert.equal(SC.formatBytes(512), '512 B')
  assert.equal(SC.formatBytes(1536), '1.5 KB')
  assert.equal(SC.formatBytes(1024 * 1024), '1 MB')
  assert.equal(SC.formatBytes(null), '—')
})

test('normalizeAnalysisResult：三种包装、figures(url/file)、obs_merged、warnings', () => {
  const direct = SC.normalizeAnalysisResult({
    ok: true,
    id: 'cell_cycle', title: '细胞周期评分', method: 'Tirosh 打分',
    summary: { bullets: ['完成打分'], metrics: { n_cells: 2400 } },
    tables: [{ id: 't1', title: '相位', columns: [{ key: 'phase', label: 'Phase', type: 'string' }], rows: [{ phase: 'G1' }] }],
    figures: [
      { id: 'a.png', title: 'UMAP', file: '/tmp/a.png', url: '/singlecell/figure?id=a.png', width: 720, height: 480 },
      { file: '/tmp/b.png' },
      { title: 'dropped' },
    ],
    obs_columns: [{ key: 'phase', kind: 'categorical' }],
    obs_merged: ['phase'],
    warnings: ['lite 实现'],
  })
  assert.equal(direct.id, 'cell_cycle')
  assert.equal(direct.method, 'Tirosh 打分')
  assert.deepEqual(direct.summary.bullets, ['完成打分'])
  assert.equal(direct.summary.metrics.n_cells, 2400)
  assert.equal(direct.tables[0].columns[0].type, 'string')
  assert.equal(direct.figures.length, 2, '没有 url/file 的图应被丢弃')
  assert.equal(direct.figures[0].url, '/singlecell/figure?id=a.png')
  assert.equal(direct.figures[1].url, '', 'file-only 图保留 file 供路径兜底')
  assert.equal(direct.figures[1].file, '/tmp/b.png')
  assert.deepEqual(direct.obsMerged, ['phase'])
  assert.deepEqual(direct.obsColumns, [{ key: 'phase', kind: 'categorical' }])
  assert.deepEqual(direct.warnings, ['lite 实现'])

  assert.equal(SC.normalizeAnalysisResult({ analysis: { id: 'wrapped' } }).id, 'wrapped')
  assert.equal(SC.normalizeAnalysisResult({ result: { id: 'wrapped2' } }).id, 'wrapped2')
  const empty = SC.normalizeAnalysisResult(null)
  assert.deepEqual([empty.summary.bullets, empty.tables, empty.figures], [[], [], []])
})

test('normalizeAnalysesPayload：groups/columns/group_keys/diagnostics', () => {
  const payload = SC.normalizeAnalysesPayload({
    dataset: 'ds1',
    groups: ['基本分析', '组间对比'],
    columns: [{ key: 'clusters' }, { key: 'n_genes', kind: 'numeric' }],
    group_keys: ['clusters', 'cell_type'],
    diagnostics: { sets: 'ok' },
    analyses: [{ id: 'cell_cycle', title: '细胞周期评分', group: '基本分析', native: true, available: true }],
  })
  assert.equal(payload.dataset, 'ds1')
  assert.deepEqual(payload.groups, ['基本分析', '组间对比'])
  assert.deepEqual(payload.columns, [{ key: 'clusters', kind: 'categorical' }, { key: 'n_genes', kind: 'numeric' }])
  assert.deepEqual(payload.groupKeys, ['clusters', 'cell_type'])
  assert.deepEqual(payload.diagnostics, { sets: 'ok' })
  assert.equal(payload.analyses.length, 1)
})

test('buildColorOptions / mergeObsColumns：固定项在前、分析新增列自动可选', () => {
  const options = SC.buildColorOptions([
    { key: 'clusters', kind: 'categorical' },
    { key: 'cell_type', kind: 'categorical' },
    { key: 'n_genes', kind: 'numeric' },
    { key: 'phase', kind: 'categorical' },
    { key: 'entropy', kind: 'numeric' },
  ])
  assert.deepEqual(options.map((option) => option.value), [
    'cluster:clusters', 'obs:cell_type', 'obs:n_genes', 'obs:pct_counts_mt',
    'obs:entropy', 'obs:phase', 'gene',
  ])
  assert.equal(options[0].labelKey, 'colorCluster')
  assert.equal(options[4].label, 'entropy · numeric')
  assert.equal(options[5].label, 'phase · categorical')
  assert.equal(options[6].labelKey, 'colorGene')

  const merged = SC.mergeObsColumns(
    [{ key: 'clusters', kind: 'categorical' }],
    [{ key: 'phase', kind: 'numeric' }, { key: 'clusters', kind: 'numeric' }],
    [{ key: 'regulon_1' }],
  )
  assert.deepEqual(merged.map((column) => column.key), ['clusters', 'phase', 'regulon_1'])
  assert.equal(merged[0].kind, 'categorical', 'preview 的 kind 优先')
  assert.equal(merged[2].kind, 'categorical', '缺 kind 时回落 categorical')
})

test('normalizeFormats：服务端清单优先、不可用标记、失败回落内置清单', () => {
  const server = SC.normalizeFormats({ formats: [
    { id: 'auto', title: '自动识别', extensions: [], available: true, note: '' },
    { id: 'rds', title: 'R 对象', extensions: ['.rds'], available: false, note: '未找到 Rscript' },
  ] })
  assert.equal(server.length, 2)
  assert.equal(server[1].available, false)
  assert.equal(server[1].note, '未找到 Rscript')

  const fallback = SC.normalizeFormats(null)
  assert.ok(fallback.length >= 8)
  assert.equal(fallback[0].id, 'auto')
  assert.ok(fallback.some((item) => item.id === 'rdata'))
  assert.ok(fallback.some((item) => item.id === 'archive'))
})

test('buildImportBody：path/url、format 同时发 data_format、member/object/transpose', () => {
  assert.deepEqual(SC.buildImportBody({ path: '/tmp/a.h5ad' }), { actor: 'user', path: '/tmp/a.h5ad' })
  assert.deepEqual(SC.buildImportBody({ url: 'https://x/y.csv.gz' }), { actor: 'user', url: 'https://x/y.csv.gz' })
  assert.deepEqual(SC.buildImportBody({ path: '/tmp/x.RData', format: 'auto' }), { actor: 'user', path: '/tmp/x.RData' })
  assert.deepEqual(
    SC.buildImportBody({ url: 'https://x/y.gz', format: 'matrix_text', member: ' a/b.csv ', object: ' seu ', transpose: 'true' }),
    {
      actor: 'user',
      url: 'https://x/y.gz',
      format: 'matrix_text',
      data_format: 'matrix_text',
      member: 'a/b.csv',
      object: 'seu',
      transpose: true,
    },
  )
  assert.equal(SC.buildImportBody({ path: '/tmp/x', transpose: 'false' }).transpose, false)
  assert.equal(SC.buildImportBody({ path: '/tmp/x', transpose: 'auto' }).transpose, undefined)
  assert.ok(SC.IMPORT_ACCEPT.includes('.rds'))
  assert.ok(SC.IMPORT_ACCEPT.includes('.RData'))
  assert.ok(SC.IMPORT_ACCEPT.includes('.zip'))
  assert.ok(SC.IMPORT_ACCEPT.includes('.h5ad'))
  assert.ok(SC.IMPORT_ACCEPT.includes('.mtx'))
})

test('normalizeImportInspect / normalizeCandidate：预检结果与候选详情', () => {
  const inspect = SC.normalizeImportInspect({
    ok: true,
    format: '10x_mtx',
    size: 2048,
    n_obs: 2400,
    n_vars: 1200,
    notes: ['自动识别'],
    candidates: [{ member: 'counts/matrix.mtx', size: 1024 }],
  })
  assert.equal(inspect.format, '10x_mtx')
  assert.equal(inspect.size, 2048)
  assert.equal(inspect.nObs, 2400)
  assert.deepEqual(inspect.notes, ['自动识别'])
  assert.equal(inspect.candidates.length, 1)
  assert.deepEqual(inspect.candidates[0].patch, { member: 'counts/matrix.mtx' })

  const candidate = SC.normalizeCandidate({ object: 'seu', class: 'Seurat', label: 'seu', n_obs: 10, n_vars: 20 }, 0)
  assert.equal(candidate.label, 'seu')
  assert.deepEqual(candidate.patch, { object: 'seu' })
  assert.match(candidate.detail, /Seurat/)
})

test('api v3：GET /formats、GET /analyses?dataset、/import/inspect、/analysis', async () => {
  const { fake, calls } = createFakeFetch([])
  const api = SC.createApi(fake)
  await api.getFormats()
  await api.getAnalyses('ds1')
  await api.inspectSource({ path: '/tmp/x.rds', format: 'rds', object: 'seu' })
  await api.importDataset({ url: 'https://x/y.csv.gz', format: 'matrix_text', member: 'm', object: 'o', transpose: 'true' })
  await api.runAnalysis('ds1', 'cell_cycle', { n_clusters: 6 })

  assert.deepEqual(calls.map((call) => call.url), [
    '/singlecell/formats',
    '/singlecell/analyses?dataset=ds1',
    '/singlecell/import/inspect',
    '/singlecell/dataset/open',
    '/singlecell/analysis',
  ])
  assert.equal(calls[0].options.method, 'GET')
  assert.equal(calls[1].options.method, 'GET')
  assert.deepEqual(JSON.parse(calls[2].options.body), {
    actor: 'user', path: '/tmp/x.rds', format: 'rds', data_format: 'rds', object: 'seu',
  })
  assert.deepEqual(JSON.parse(calls[3].options.body), {
    actor: 'user',
    url: 'https://x/y.csv.gz',
    format: 'matrix_text',
    data_format: 'matrix_text',
    member: 'm',
    object: 'o',
    transpose: true,
  })
  assert.deepEqual(JSON.parse(calls[4].options.body), {
    dataset: 'ds1', id: 'cell_cycle', params: { n_clusters: 6 }, actor: 'user',
  })
})

test('api ambiguous：错误对象带 code/candidates/status，供候选列表重试', async () => {
  const { fake } = createFakeFetch([{
    body: {
      ok: false,
      code: 'ambiguous',
      error: '检测到多个候选对象',
      candidates: [{ object: 'seu', class: 'Seurat', n_obs: 10, n_vars: 20 }],
    },
    status: 400,
  }])
  const api = SC.createApi(fake)
  await assert.rejects(
    () => api.importDataset({ path: '/tmp/x.RData' }),
    (error) => {
      assert.equal(error.code, 'ambiguous')
      assert.equal(error.status, 400)
      assert.equal(error.candidates.length, 1)
      assert.equal(error.candidates[0].object, 'seu')
      assert.match(error.message, /多个候选/)
      return true
    },
  )
})

// ---------------------------------------------------------------------------
// 8.6 v4：样本 & 分组（GEO 每个 GSM 一个文件 → 合并 + 改名 + 大组）
// ---------------------------------------------------------------------------

test('parsePathLines：多行/逗号分隔、去空行与重复', () => {
  assert.deepEqual(SC.parsePathLines(' /a/x.csv.gz \n/b/y.csv.gz\n\n/a/x.csv.gz, /c '), [
    '/a/x.csv.gz', '/b/y.csv.gz', '/c',
  ])
  assert.deepEqual(SC.parsePathLines('/a，/b'), ['/a', '/b'])
  assert.deepEqual(SC.parsePathLines(''), [])
  assert.deepEqual(SC.parsePathLines(null), [])
})

test('normalizeSamples / normalizeObsPayload：兼容 Host 真实返回', () => {
  const samples = SC.normalizeSamples([
    { id: 'young1', name: '年轻1', group: 'young', n_cells: 300 },
    { id: 'old1', name: '老年1', group: 'old' },
    { name: '没有 id 的样本' },
  ])
  assert.equal(samples.length, 2)
  assert.equal(samples[0].nCells, 300)
  assert.equal(samples[0].name, '年轻1')
  assert.equal(samples[1].nCells, 0)

  const obs = SC.normalizeObsPayload({
    ok: true,
    dataset: 'ds1',
    columns: [
      { key: 'sample', kind: 'categorical', n_unique: 3, values: ['young1', 'old1'] },
      { key: 'n_genes', kind: 'numeric', n_unique: 100 },
    ],
    samples: [{ id: 'young1', name: '年轻1', group: 'young', n_cells: 300 }],
    groups: ['old', 'young'],
    group_column: 'group',
    sample_column: 'sample',
  })
  assert.equal(obs.dataset, 'ds1')
  assert.equal(obs.columns.length, 2)
  assert.deepEqual(obs.columns[0].values, ['young1', 'old1'])
  assert.equal(obs.columns[1].kind, 'numeric')
  assert.equal(obs.samples[0].group, 'young')
  assert.deepEqual(obs.groups, ['old', 'young'])
  assert.equal(obs.groupColumn, 'group')
  assert.equal(obs.sampleColumn, 'sample')
  assert.deepEqual(SC.normalizeObsPayload(null).groups, [])
})

test('parseSampleRows / groupOptions：样本表与分组下拉', () => {
  const payload = {
    samples: [
      { id: 'young1', name: '年轻1', group: 'young', n_cells: 300 },
      { id: 'old1', name: '老年1', group: 'old', n_cells: 300 },
    ],
    groups: ['old', 'young'],
  }
  const rows = SC.parseSampleRows(payload)
  assert.deepEqual(rows.map((row) => row.id), ['young1', 'old1'])
  assert.equal(rows[0].name, '年轻1')
  assert.equal(rows[0].dirty, false)

  const options = SC.groupOptions({ ...payload, samples: payload.samples.concat([{ id: 'm', name: 'm', group: 'mid' }]) })
  assert.deepEqual(options.map((option) => option.value), ['', 'mid', 'old', 'young'])
})

test('buildSamplesUpdate / applySampleGroupMapping：全量提交 + 预览映射应用', () => {
  const rows = [
    { id: 'young1', name: '年轻1', group: 'young', nCells: 300 },
    { id: 'old1', name: '', group: '' },
  ]
  const body = SC.buildSamplesUpdate(rows, { groups: ['young'], sampleColumn: 'sample', groupColumn: 'group' })
  assert.deepEqual(body.samples, [
    { id: 'young1', name: '年轻1', group: 'young' },
    { id: 'old1', group: '' },
  ])
  assert.deepEqual(body.groups, ['young'])
  assert.equal(body.sample_column, 'sample')
  assert.equal(body.group_column, 'group')

  const mapped = SC.applySampleGroupMapping(rows, { young1: 'old', old1: 'young' })
  assert.equal(mapped[0].group, 'old')
  assert.equal(mapped[0].dirty, true)
  assert.equal(mapped[1].group, 'young')
  assert.equal(mapped[1].name, '')
  assert.equal(rows[0].group, 'young', '原数组不应被修改')
})

test('auto_group / import-map 结果归一化与预览', () => {
  const auto = SC.normalizeAutoGroupResult({ mapping: { young1: 'young', old1: 'old', bad: '' }, unmatched: ['x'] })
  assert.deepEqual(auto.mapping, { young1: 'young', old1: 'old' })
  assert.deepEqual(auto.unmatched, ['x'])
  const preview = SC.parseAutoGroupMapping({ mapping: { old1: 'old', young1: 'young' } })
  assert.deepEqual(preview.rows, [
    { id: 'old1', group: 'old' },
    { id: 'young1', group: 'young' },
  ])

  const mapped = SC.normalizeImportMapResult({
    updates: [{ id: 'young1', group: 'young' }, { id: 'old1', group: 'old' }],
    matched: 2,
    total: 3,
  })
  assert.deepEqual(mapped.updates, { young1: 'young', old1: 'old' })
  assert.equal(mapped.matched, 2)
  assert.equal(mapped.total, 3)
})

test('obs join 结果与匹配率文案', () => {
  const result = SC.normalizeJoinObsResult({ joined_columns: ['age', 'sex'], matched: 850, total: 900, match_rate: 850 / 900 })
  assert.deepEqual(result.joinedColumns, ['age', 'sex'])
  assert.equal(result.matched, 850)
  assert.equal(SC.formatMatchRate(result), '已匹配 850/900（94.4%）')
  assert.equal(SC.formatMatchRate({ matched: 3, total: 3 }), '已匹配 3/3（100%）')
  assert.equal(SC.formatMatchRate({}), '')
})

test('describeGroupState：当前分组列与组列表', () => {
  assert.equal(SC.describeGroupState(['old', 'young'], 'group'), '当前分组列：group（2 组：old、young）')
  assert.equal(SC.describeGroupState([], 'group'), '当前分组列：group（未分组）')
})

test('buildImportBody v4：paths 多文件与 merge_dir 目录合并', () => {
  assert.deepEqual(SC.buildImportBody({ paths: ['/a/x.csv.gz', '/a/y.csv.gz'] }), {
    actor: 'user',
    paths: ['/a/x.csv.gz', '/a/y.csv.gz'],
  })
  assert.deepEqual(SC.buildImportBody({ paths: ['/a/x.csv.gz'] }), { actor: 'user', path: '/a/x.csv.gz' })
  assert.deepEqual(SC.buildImportBody({ path: '/data/geo', mergeDir: true }), {
    actor: 'user',
    path: '/data/geo',
    merge_dir: true,
  })
})

test('paramDefaults v4：column 参数在分组列存在时默认落到 group_column', () => {
  const params = [
    { key: 'group', type: 'column', default: '' },
    { key: 'cell_type', type: 'column', default: 'cell_type' },
    { key: 'n', type: 'number', default: 5 },
  ]
  const defaults = SC.paramDefaults(params, [{ key: 'clusters', kind: 'categorical' }, { key: 'group', kind: 'categorical' }], 'group')
  assert.equal(defaults.group, 'group')
  assert.equal(defaults.cell_type, 'cell_type', '不存在的列保持原默认值')
  assert.equal(defaults.n, 5)
  // 旧签名（无 obs/分组列）行为不变
  assert.equal(SC.paramDefaults(params).group, '')
})

test('normalizeDataset v4：samples / groups / group_column 进 state 数据', () => {
  const dataset = SC.normalizeDataset({
    id: 'ds9', name: 'merged-3-samples', n_obs: 900, n_vars: 1200, steps: [],
    samples: [{ id: 'young1', name: '年轻1', group: 'young', n_cells: 300 }],
    groups: ['old', 'young'],
    group_column: 'group',
    sample_column: 'sample',
  })
  assert.equal(dataset.samples.length, 1)
  assert.equal(dataset.samples[0].nCells, 300)
  assert.deepEqual(dataset.groups, ['old', 'young'])
  assert.equal(dataset.groupColumn, 'group')
  assert.equal(dataset.sampleColumn, 'sample')
})

test('api v4：/obs、/samples/update、/samples/auto-group、/samples/import-map、/obs/join、/obs/set-group、/obs/derive', async () => {
  const { fake, calls } = createFakeFetch([])
  const api = SC.createApi(fake)
  await api.getObs('ds1')
  await api.samplesUpdate('ds1', {
    samples: [{ id: 'young1', name: '年轻1', group: 'young' }],
    groups: ['young', 'old'],
    sample_column: 'sample',
    group_column: 'group',
  })
  await api.samplesAutoGroup('ds1', '(GC|LM)')
  await api.samplesImportMap('ds1', '/tmp/map.csv')
  await api.obsJoin('ds1', '/tmp/meta.csv', 'cell')
  await api.obsSetGroup('ds1', 'group')
  await api.obsDerive('ds1', { name: 'batch', source: 'cell_name', pattern: '^(\\w+)_' })

  assert.deepEqual(calls.map((call) => call.url), [
    '/singlecell/obs?dataset=ds1',
    '/singlecell/samples/update',
    '/singlecell/samples/auto-group',
    '/singlecell/samples/import-map',
    '/singlecell/obs/join',
    '/singlecell/obs/set-group',
    '/singlecell/obs/derive',
  ])
  assert.equal(calls[0].options.method, 'GET')
  assert.deepEqual(JSON.parse(calls[1].options.body), {
    dataset: 'ds1',
    actor: 'user',
    samples: [{ id: 'young1', name: '年轻1', group: 'young' }],
    groups: ['young', 'old'],
    sample_column: 'sample',
    group_column: 'group',
  })
  assert.deepEqual(JSON.parse(calls[2].options.body), { dataset: 'ds1', pattern: '(GC|LM)', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[3].options.body), { dataset: 'ds1', path: '/tmp/map.csv', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[4].options.body), { dataset: 'ds1', path: '/tmp/meta.csv', key: 'cell', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[5].options.body), { dataset: 'ds1', column: 'group', actor: 'user' })
  assert.deepEqual(JSON.parse(calls[6].options.body), { dataset: 'ds1', actor: 'user', name: 'batch', source: 'cell_name', pattern: '^(\\w+)_' })
})

test('静态：v4 样本 & 分组 UI 与接口接线', () => {
  assert.match(CLIENT_SOURCE, /function SamplesPanel\(props\)/)
  assert.match(CLIENT_SOURCE, /samplesPanel: samplesPanelNode/)
  assert.match(CLIENT_SOURCE, /className: 'sc-samples'/)
  assert.match(CLIENT_SOURCE, /parsePathLines\(path\)/)
  assert.match(CLIENT_SOURCE, /mergeDir/)
  assert.match(CLIENT_SOURCE, /importMergedPrefix/)
  for (const route of ["'/samples/update'", "'/samples/auto-group'", "'/samples/import-map'", "'/obs/join'", "'/obs/set-group'", "'/obs/derive'", "'/obs' + query"]) {
    assert.ok(CLIENT_SOURCE.includes(route), `缺少 v4 接口接线：${route}`)
  }
  assert.match(CLIENT_SOURCE, /\.sc-samples-table/)
  assert.match(CLIENT_SOURCE, /\.sc-path-lines/)
})

// ---------------------------------------------------------------------------
// 8.7 v5：数据根目录扫描 → 样本表 → 确认导入（GEO 输入逻辑）
// ---------------------------------------------------------------------------

test('normalizeScanResult：扫描结果归一化（root/samples/errors）', () => {
  const scan = SC.normalizeScanResult({
    ok: true,
    root: '/data/GSE246662_RAW',
    samples: [
      { id: 's1', name: 'GSM1_GC3', raw_name: 'GSM1_GC3', format: 'expr_matrix', format_label: '表达矩阵', path: '/data/a/GSM1_GC3.csv.gz', detail: '逗号分隔', trio_prefix: null },
      { path: '/data/b', format: '10x_mtx' },
    ],
    errors: ['子目录 x 无法识别'],
  })
  assert.equal(scan.root, '/data/GSE246662_RAW')
  assert.equal(scan.samples.length, 2)
  assert.equal(scan.samples[0].formatLabel, '表达矩阵')
  assert.equal(scan.samples[0].trioPrefix, '')
  assert.equal(scan.samples[1].name, 'b', '缺 name 时用路径 basename')
  assert.deepEqual(scan.errors, ['子目录 x 无法识别'])
})

test('createScanRows：默认全选、未分组、key 唯一', () => {
  const rows = SC.createScanRows({
    root: '/r',
    samples: [
      { name: 'S1', raw_name: 'S1', format: 'h5ad', path: '/r/a.h5ad' },
      { name: 'S1', raw_name: 'S1', format: 'h5ad', path: '/r/a.h5ad' },
    ],
  })
  assert.equal(rows.length, 2)
  assert.equal(rows[0].selected, true)
  assert.equal(rows[0].group, '')
  assert.notEqual(rows[0].key, rows[1].key)
  assert.deepEqual(SC.createScanRows(null), [])
})

test('findDuplicateSampleNames / validateSampleSelection：改名去重与开始分析校验', () => {
  assert.deepEqual(SC.findDuplicateSampleNames([{ name: 'A' }, { name: 'a ' }, { name: 'B' }]), ['A'])
  const rows = [
    { key: '1', selected: true, name: 'PT1', group: 'young', rawName: 'GSM1', path: '/a' },
    { key: '2', selected: true, name: 'PT1', group: 'old', rawName: 'GSM2', path: '/b' },
    { key: '3', selected: false, name: '', group: '', rawName: 'GSM3', path: '/c' },
  ]
  const invalid = SC.validateSampleSelection(rows)
  assert.equal(invalid.ok, false)
  assert.ok(invalid.errors.some((item) => item.includes('样本名重复')))
  const noGroup = SC.validateSampleSelection([
    { key: '1', selected: true, name: 'PT1', group: '', rawName: 'GSM1', path: '/a' },
    { key: '2', selected: true, name: 'PT2', group: 'old', rawName: 'GSM2', path: '/b' },
  ])
  assert.equal(noGroup.ok, false)
  assert.ok(noGroup.errors.join('').includes('没有分组'))
  assert.equal(SC.validateSampleSelection([{ key: '1', selected: false, name: 'x', group: 'y' }]).ok, false)
  const valid = SC.validateSampleSelection([
    { key: '1', selected: true, name: 'PT1', group: 'young', rawName: 'GSM1', path: '/a' },
    { key: '2', selected: true, name: 'PT2', group: 'old', rawName: 'GSM2', path: '/b' },
  ])
  assert.equal(valid.ok, true)
  assert.equal(valid.selectedCount, 2)
})

test('buildSamplesOpenPayload：只送勾选样本，带 name/group/trio_prefix/format', () => {
  const payload = SC.buildSamplesOpenPayload([
    { key: '1', selected: true, name: 'PT1', group: 'young', rawName: 'GSM1', path: '/a/GSM1.csv.gz', format: 'expr_matrix' },
    { key: '2', selected: true, name: '  ', group: 'old', rawName: 'GSM2_RAW', path: '/a/GSM2', trioPrefix: 'GSM2' },
    { key: '3', selected: false, name: 'PT3', group: 'old', rawName: 'GSM3', path: '/a/GSM3' },
  ])
  assert.deepEqual(payload.samples, [
    { path: '/a/GSM1.csv.gz', name: 'PT1', group: 'young', format: 'expr_matrix' },
    { path: '/a/GSM2', name: 'GSM2_RAW', group: 'old', trio_prefix: 'GSM2' },
  ])
  assert.equal(payload.name, undefined)
  assert.equal(SC.buildSamplesOpenPayload([], { name: 'my-dataset' }).name, 'my-dataset')
})

test('setScanRowsSelected / applyGroupToRows：全选/清空与批量分组', () => {
  const rows = [
    { key: '1', selected: true, group: '' },
    { key: '2', selected: false, group: 'old' },
  ]
  assert.deepEqual(SC.setScanRowsSelected(rows, false).map((row) => row.selected), [false, false])
  assert.deepEqual(SC.setScanRowsSelected(rows, true).map((row) => row.selected), [true, true])
  const grouped = SC.applyGroupToRows(rows, 'young', true)
  assert.equal(grouped[0].group, 'young')
  assert.equal(grouped[1].group, 'old', '未勾选行不受影响')
  assert.equal(SC.applyGroupToRows(rows, 'x', false)[1].group, 'x')
  assert.equal(rows[0].group, '', '原数组不变')
})

test('normalizeGeoResult / matchGeoToSamples：GSM 说明按文件名匹配到样本行', () => {
  const geo = SC.normalizeGeoResult({
    gse: 'GSE246662',
    title: 'Test GEO',
    samples: [
      { gsm: 'GSM7874174', title: 'GC3 sample' },
      { gsm: 'GSM7874175', title: 'GC4 sample' },
    ],
  })
  assert.equal(geo.title, 'Test GEO')
  assert.equal(geo.samples.length, 2)
  const matched = SC.matchGeoToSamples([
    { key: '1', name: 'GSM7874174_GC3', rawName: 'GSM7874174_GC3', trioPrefix: '', path: '/d/GSM7874174_GC3.csv.gz', geoTitle: '' },
    { key: '2', name: 'other', rawName: 'other', trioPrefix: 'GSM7874175_x', path: '/d/other', geoTitle: '' },
    { key: '3', name: 'none', rawName: 'none', trioPrefix: '', path: '/d/none', geoTitle: '' },
  ], { gse: 'GSE246662', samples: [{ gsm: 'GSM7874174', title: 'GC3 sample' }, { gsm: 'GSM7874175', title: 'GC4 sample' }] })
  assert.equal(matched.rows[0].geoTitle, 'GC3 sample')
  assert.equal(matched.rows[1].geoTitle, 'GC4 sample')
  assert.equal(matched.rows[2].geoTitle, '')
  assert.equal(matched.matched, 2)
  assert.equal(matched.unmatched.length, 0)
  const partial = SC.matchGeoToSamples([{ key: '1', rawName: 'nothing', path: '/x' }], {
    samples: [{ gsm: 'GSM1', title: 'a' }, { gsm: 'GSM2', title: 'b' }],
  })
  assert.equal(partial.unmatched.length, 2)
})

test('normalizeImportState / formatImportStateHint：已跑过一半的数据提示', () => {
  const state = SC.normalizeImportState({
    steps: ['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'embed', 'bogus'],
    current: 'embed',
    details: {},
    mapped: { cluster_column: 'seurat_clusters', celltype_column: 'cell_type', n_clusters: 8, n_variable_features: 2000, has_pca: true, has_embed: true },
  })
  assert.deepEqual(state.steps, ['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'embed'])
  assert.equal(state.mapped.clusterColumn, 'seurat_clusters')
  assert.equal(state.mapped.hasEmbed, true)
  const hint = SC.formatImportStateHint(state, ['已有细胞类型列：cell_type'])
  assert.ok(hint.includes('检测到已有进度：已完成 8 步'))
  assert.ok(hint.includes('聚类列 seurat_clusters（8 簇）'))
  assert.ok(hint.includes('细胞类型列 cell_type'))
  assert.ok(hint.includes('已有二维嵌入'))
  assert.ok(hint.includes('可继续往下跑'))
  assert.equal(SC.formatImportStateHint(null, []), '')
})

test('buildImportBody v5：samples 载荷（确认样本与分组后 open）', () => {
  const body = SC.buildImportBody({
    samples: [
      { path: '/a/GSM1.csv.gz', name: 'PT1', group: 'young', format: 'expr_matrix' },
      { path: '/a/GSM2.csv.gz', name: 'PT2', group: 'old', trio_prefix: 'GSM2' },
      { name: '缺 path' },
    ],
    name: 'geo-merged',
  })
  assert.deepEqual(JSON.parse(JSON.stringify(body)), {
    actor: 'user',
    samples: [
      { path: '/a/GSM1.csv.gz', name: 'PT1', group: 'young', format: 'expr_matrix' },
      { path: '/a/GSM2.csv.gz', name: 'PT2', group: 'old', trio_prefix: 'GSM2' },
    ],
    name: 'geo-merged',
  })
})

test('normalizePreview v5：带 notes 与 import_state', () => {
  const preview = SC.normalizePreview({
    genes: ['CD3D'],
    obs_columns: [{ key: 'clusters', kind: 'categorical' }],
    n_obs: 100,
    n_vars: 200,
    steps: ['qc', 'pca'],
    notes: ['检测到已有进度'],
    import_state: { steps: ['qc', 'pca'], current: 'pca', mapped: { has_pca: true } },
  })
  assert.deepEqual(preview.notes, ['检测到已有进度'])
  assert.deepEqual(preview.importState.steps, ['qc', 'pca'])
  assert.equal(preview.importState.mapped.hasPca, true)
})

test('api v5：POST /scan 与 POST /geo', async () => {
  const { fake, calls } = createFakeFetch([])
  const api = SC.createApi(fake)
  await api.scanDataRoot('/data/GSE1_RAW', true)
  await api.queryGeo('gse246662')
  assert.deepEqual(calls.map((call) => call.url), ['/singlecell/scan', '/singlecell/geo'])
  assert.equal(calls[0].options.method, 'POST')
  assert.deepEqual(JSON.parse(calls[0].options.body), { path: '/data/GSE1_RAW', recursive: true, actor: 'user' })
  assert.deepEqual(JSON.parse(calls[1].options.body), { gse: 'GSE246662', actor: 'user' })
})

test('静态：v5 扫描样本表 UI 与 import_state 提示接线', () => {
  assert.match(CLIENT_SOURCE, /function normalizeScanResult\(payload\)/)
  assert.match(CLIENT_SOURCE, /function validateSampleSelection\(rows\)/)
  assert.match(CLIENT_SOURCE, /function buildSamplesOpenPayload\(rows, options\)/)
  assert.match(CLIENT_SOURCE, /function formatImportStateHint\(importState, notes\)/)
  assert.match(CLIENT_SOURCE, /scanDataRoot\(root, scanRecursive !== false\)/)
  assert.match(CLIENT_SOURCE, /queryGeo\(gse\)/)
  assert.match(CLIENT_SOURCE, /onConfirmSamples: handleConfirmSamples/)
  assert.match(CLIENT_SOURCE, /importStateHintText/)
  assert.match(CLIENT_SOURCE, /scanAdvancedToggle/)
  for (const route of ["'/scan'", "'/geo'"]) {
    assert.ok(CLIENT_SOURCE.includes(route), `缺少 v5 接口接线：${route}`)
  }
})

// ---------------------------------------------------------------------------
// 8.8 v6：降维聚类控制卡 / 手动注释 / 分组维护 / Agent 意图解析
// ---------------------------------------------------------------------------

test('normalizeDimredSettings / clampResolution：默认值与范围', () => {
  const defaults = SC.normalizeDimredSettings(null)
  assert.equal(defaults.engine, 'seurat')
  assert.equal(defaults.harmony, true)
  assert.equal(defaults.resolution, 0.2)
  assert.equal(defaults.nPcs, 30)
  assert.equal(defaults.neighbors, 30)
  assert.equal(defaults.nTopGenes, 3000)

  assert.equal(SC.clampResolution(2), 1.5)
  assert.equal(SC.clampResolution(0), 0.05)
  assert.equal(SC.clampResolution(0.333), 0.33)
  assert.equal(SC.clampResolution('abc'), 0.2)
  const custom = SC.normalizeDimredSettings({ engine: 'builtin', harmony: false, resolution: 0.8, nPcs: 50, neighbors: 20, nTopGenes: 2000 })
  assert.equal(custom.engine, 'builtin')
  assert.equal(custom.harmony, false)
  assert.equal(custom.resolution, 0.8)
  assert.equal(custom.nPcs, 50)
  assert.equal(custom.neighbors, 20)
  assert.equal(custom.nTopGenes, 2000)
})

test('buildReclusterParams / buildReclusterEmbedParams：Seurat+Harmony+batch_key payload', () => {
  const settings = SC.normalizeDimredSettings({
    engine: 'seurat', harmony: true, batchKey: '', resolution: 0.3, nPcs: 30, neighbors: 30,
  })
  const cluster = SC.buildReclusterParams(settings, 'sample')
  assert.equal(cluster.engine, 'seurat')
  assert.equal(cluster.use_harmony, true)
  assert.equal(cluster.resolution, 0.3)
  assert.equal(cluster.n_pcs, 30)
  assert.equal(cluster.n_neighbors, 30)
  assert.equal(cluster.batch_key, 'sample')
  assert.equal(cluster.method, 'louvain')

  const explicit = SC.buildReclusterParams(SC.normalizeDimredSettings({ batchKey: 'orig.ident' }), 'sample')
  assert.equal(explicit.batch_key, 'orig.ident')

  const builtin = SC.buildReclusterParams(SC.normalizeDimredSettings({ engine: 'builtin' }), 'sample')
  assert.equal(builtin.engine, undefined)
  assert.equal(builtin.use_harmony, undefined)
  assert.equal(builtin.batch_key, undefined)
  assert.equal(builtin.resolution, 0.2)

  const embed = SC.buildReclusterEmbedParams(settings, 'sample')
  assert.equal(embed.layout, 'umap')
  assert.equal(embed.engine, 'seurat')
  assert.equal(embed.use_harmony, true)
  assert.equal(embed.batch_key, 'sample')
  const builtinEmbed = SC.buildReclusterEmbedParams(SC.normalizeDimredSettings({ engine: 'builtin' }))
  assert.equal(builtinEmbed.engine, undefined)
  assert.equal(builtinEmbed.layout, 'umap')
  assert.equal(SC.buildDimredRecipeRequest({}).id, 'seurat6')
  assert.equal(SC.DIMRED_RECIPE_ID, 'seurat6')
})

test('pickDimredWarnings：从 dataset warnings 与运行 bullets 里挑提示', () => {
  const warnings = SC.pickDimredWarnings(
    ['聚类数 35 偏多（参考 20–30）：建议把 resolution 调低（当前 0.80）', '其它无关警告'],
    [{ summary: { bullets: ['完成 Seurat 6.0 聚类', 'Harmony 整合：3 个批次'] } }],
  )
  assert.ok(warnings.some((item) => item.indexOf('聚类数 35 偏多') >= 0))
  assert.ok(warnings.some((item) => item.indexOf('Harmony') >= 0))
  assert.ok(!warnings.includes('其它无关警告'))
})

test('buildClusterAnnotationRows / parseAnnotationCsv：手动注释映射', () => {
  const rows = SC.buildClusterAnnotationRows([
    { id: '0', label: '0', nCells: 312, genes: [{ gene: 'CD3D' }, { gene: 'CD3E' }, { gene: 'IL7R' }, { gene: 'X' }] },
    { id: '1', label: '1', nCells: 200, genes: [] },
  ], { 0: 'T cell' })
  assert.equal(rows.length, 2)
  assert.deepEqual(rows[0].topGenes, ['CD3D', 'CD3E', 'IL7R'])
  assert.equal(rows[0].cellType, 'T cell')
  assert.equal(rows[1].cellType, '')

  const parsed = SC.parseAnnotationCsv('cluster,cell_type\n0,T cell\n1,B cell\n2\n')
  assert.deepEqual(parsed.mapping, { 0: 'T cell', 1: 'B cell' })
  assert.equal(parsed.rows.length, 2)
  assert.equal(parsed.skipped, 1)
  assert.deepEqual(SC.parseAnnotationCsv('').mapping, {})
})

test('isGroupUsedBySamples / removeGroupFromList / clearSampleGroupValue / clearSampleRowGroup', () => {
  const rows = [
    { id: 's1', name: 'a', group: 'young' },
    { id: 's2', name: 'b', group: 'old' },
    { id: 's3', name: 'c', group: '' },
  ]
  assert.equal(SC.isGroupUsedBySamples(rows, 'young'), true)
  assert.equal(SC.isGroupUsedBySamples(rows, 'mid'), false)
  assert.deepEqual(SC.removeGroupFromList(['young', 'old', 'mid'], 'old'), ['young', 'mid'])
  assert.equal(SC.clearSampleGroupValue('__clear__'), '')
  assert.equal(SC.clearSampleGroupValue('young'), 'young')
  const cleared = SC.clearSampleRowGroup(rows, 's1')
  assert.equal(cleared[0].group, '')
  assert.equal(cleared[1].group, 'old')
  assert.equal(rows[0].group, 'young', '原数组不变')
})

test('extractResolution / extractGene：自然语言里的数值与基因', () => {
  assert.equal(SC.extractResolution('resolution 调到 0.3'), 0.3)
  assert.equal(SC.extractResolution('把分辨率改成 0.8'), 0.8)
  assert.equal(SC.extractResolution('res=1.2'), 1.2)
  assert.equal(SC.extractResolution('加 harmony'), null)
  assert.equal(SC.extractGene('画 CD3D 的 featureplot'), 'CD3D')
  assert.equal(SC.extractGene('看看 IL7R 表达'), 'IL7R')
  assert.equal(SC.extractGene('画 UMAP 图'), '')
  assert.equal(SC.extractGene('导出 CSV'), '')
})

test('parseAgentIntent：harmony/resolution/marker/featureplot/注释/组间/自检/导出', () => {
  const ctx = { steps: ['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'embed'], hasMarkers: true }
  const harmony = SC.parseAgentIntent('加 harmony 整合一下', ctx)
  assert.equal(harmony.fallback, false)
  assert.ok(harmony.suggestions.some((item) => item.action.type === 'recluster' && item.action.params.use_harmony === true))

  const resolution = SC.parseAgentIntent('resolution 调到 0.3', ctx)
  assert.ok(resolution.suggestions.some((item) => item.action.type === 'recluster' && item.action.params.resolution === 0.3))

  const markers = SC.parseAgentIntent('帮我找 marker', ctx)
  assert.ok(markers.suggestions.some((item) => item.action.type === 'markers'))

  const feature = SC.parseAgentIntent('把 CD3D 画到 UMAP 上', ctx)
  assert.ok(feature.suggestions.some((item) => item.action.type === 'featureplot' && item.action.gene === 'CD3D'))

  const annotate = SC.parseAgentIntent('把 cluster 0 标成 T cell', ctx)
  assert.ok(annotate.suggestions.some((item) => item.action.type === 'annotate'))

  const de = SC.parseAgentIntent('跑一下组间差异基因', ctx)
  assert.ok(de.suggestions.some((item) => item.action.type === 'analysis' && item.action.id === 'de_group'))

  const audit = SC.parseAgentIntent('做个自检', ctx)
  assert.ok(audit.suggestions.some((item) => item.action.type === 'audit'))

  const report = SC.parseAgentIntent('导出 report', ctx)
  assert.ok(report.suggestions.some((item) => item.action.type === 'export' && item.action.what === 'report'))

  const unknown = SC.parseAgentIntent('今天天气怎么样', { steps: [] })
  assert.equal(unknown.fallback, true)
  assert.equal(unknown.suggestions.length, 0)

  // 2–4 条：至少补一条上下文建议
  const single = SC.parseAgentIntent('加 harmony', { steps: [] })
  assert.ok(single.suggestions.length >= 2)
  assert.ok(SC.parseAgentIntent('加 harmony', ctx).suggestions.length <= 4)
})

test('agentActionInstruction：建议动作 → 中文指令', () => {
  assert.match(SC.agentActionInstruction({ type: 'recipe', id: 'seurat6' }), /一键降维聚类/)
  const recluster = SC.agentActionInstruction({ type: 'recluster', params: { engine: 'seurat', use_harmony: true, resolution: 0.3 } })
  assert.match(recluster, /只重跑聚类/)
  assert.match(recluster, /Seurat 引擎/)
  assert.match(recluster, /Harmony/)
  assert.match(recluster, /resolution=0\.3/)
  assert.match(SC.agentActionInstruction({ type: 'markers' }), /标记基因/)
  assert.match(SC.agentActionInstruction({ type: 'featureplot', gene: 'CD3D' }), /CD3D/)
  assert.match(SC.agentActionInstruction({ type: 'analysis', id: 'de_group' }), /de_group/)
  assert.match(SC.agentActionInstruction({ type: 'export', what: 'methods' }), /methods/)
  assert.match(SC.agentActionInstruction(null), /继续/)
})

test('静态：v6 降维卡 / 注释工作台 / 分组删除 / Agent 输入接线', () => {
  assert.match(CLIENT_SOURCE, /function DimredCard\(props\)/)
  assert.match(CLIENT_SOURCE, /function AnnotationWorkbench\(props\)/)
  assert.match(CLIENT_SOURCE, /onRunPreset: handleStep2Preset/)
  assert.match(CLIENT_SOURCE, /onRecluster: handleStep2Recluster/)
  assert.match(CLIENT_SOURCE, /buildReclusterParams\(settings, batchKey\)/)
  assert.match(CLIENT_SOURCE, /buildReclusterEmbedParams\(settings, batchKey\)/)
  assert.match(CLIENT_SOURCE, /onPlotMarker: handlePlotMarker/)
  assert.match(CLIENT_SOURCE, /onFeaturePlot: handleFeaturePlot/)
  assert.match(CLIENT_SOURCE, /parseAnnotationCsv\(annotCsvText\)/)
  assert.match(CLIENT_SOURCE, /onDeleteGroup: handleDeleteSampleGroup/)
  assert.match(CLIENT_SOURCE, /onRemoveSample: handleRemoveSample/)
  assert.match(CLIENT_SOURCE, /onAgentRun: handleAgentRun/)
  assert.match(CLIENT_SOURCE, /agentAskPlaceholder/)
  assert.match(CLIENT_SOURCE, /parseAgentIntent\(text/)
  assert.match(CLIENT_SOURCE, /DIMRED_RECIPE_ID/)
  assert.match(CLIENT_SOURCE, /sc-dimred/)
})

test('分组清空/删除 与 buildSamplesUpdate 组合：删空组后 update 载荷不含该组', () => {
  const rows = [
    { id: 's1', name: 'PT1', group: 'young', nCells: 100 },
    { id: 's2', name: 'PT2', group: '', nCells: 100 },
  ]
  assert.equal(SC.isGroupUsedBySamples(rows, 'old'), false)
  const nextGroups = SC.removeGroupFromList(['young', 'old'], 'old')
  const body = SC.buildSamplesUpdate(rows, { groups: nextGroups, sampleColumn: 'sample', groupColumn: 'group' })
  assert.deepEqual(body.groups, ['young'])
  assert.equal(body.samples[1].group, '')
  const cleared = SC.applyGroupToRows(rows, '', false)
  assert.deepEqual(cleared.map((row) => row.group), ['', ''])
})

test('parseAgentIntent：比例/通讯 与 未知句子的 fallback', () => {
  const proportion = SC.parseAgentIntent('看组间细胞比例差异', { steps: ['cluster'] })
  assert.ok(proportion.suggestions.some((item) => item.action.type === 'analysis' && item.action.id === 'proportion_test'))
  const lr = SC.parseAgentIntent('做细胞通讯', { steps: ['cluster'] })
  assert.ok(lr.suggestions.some((item) => item.action.type === 'analysis' && item.action.id === 'lr_communication'))
  const methods = SC.parseAgentIntent('导出 methods', { steps: ['cluster'] })
  assert.ok(methods.suggestions.some((item) => item.action.type === 'export' && item.action.what === 'methods'))
  assert.deepEqual(SC.parseAgentIntent('', {}), { suggestions: [], fallback: true })
})

// ---------------------------------------------------------------------------
// 8.9 v7：新 Shell（步骤导航 / Step2 参数 / 三套结果 / 日志 / 画廊 / 最近数据）
// ---------------------------------------------------------------------------

test('STEP_NAV：4 个步骤的顺序与默认步骤推导', () => {
  assert.deepEqual(SC.STEP_NAV.map((step) => step.id), ['read', 'cluster', 'annotate', 'results'])
  assert.deepEqual(SC.STEP_NAV_IDS, ['read', 'cluster', 'annotate', 'results'])
  assert.equal(SC.normalizeStepNavId('cluster'), 'cluster')
  assert.equal(SC.normalizeStepNavId('nope'), 'read')
  assert.equal(SC.defaultStepForDataset(null), 'read')
  assert.equal(SC.defaultStepForDataset({ id: 'ds', steps: [] }), 'read')
  assert.equal(SC.defaultStepForDataset({ id: 'ds', steps: ['qc', 'filter'] }), 'cluster')
  assert.equal(SC.defaultStepForDataset({ id: 'ds', steps: ['qc', 'filter', 'cluster'] }), 'annotate')
  assert.equal(SC.defaultStepForDataset({ id: 'ds', steps: ['qc', 'filter', 'cluster', 'annotate'] }), 'results')
})

test('computeStepNavStatuses：未开始 / 进行中 / 已完成 / warning', () => {
  const empty = SC.computeStepNavStatuses(null)
  assert.equal(empty.read.status, 'empty')
  assert.equal(empty.cluster.status, 'empty')
  assert.equal(empty.annotate.status, 'empty')
  assert.equal(empty.results.status, 'empty')

  const partial = SC.computeStepNavStatuses({ id: 'ds', steps: ['qc'], warnings: [] })
  assert.equal(partial.read.status, 'partial')
  assert.equal(partial.cluster.status, 'empty')

  const done = SC.computeStepNavStatuses({
    id: 'ds',
    steps: ['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'embed', 'markers', 'annotate'],
    warnings: [],
    analysesRun: [{ id: 'de_group' }],
  })
  assert.equal(done.read.status, 'done')
  assert.equal(done.cluster.status, 'done')
  assert.equal(done.annotate.status, 'done')
  assert.equal(done.results.status, 'done')
  assert.equal(done.cluster.warn, false)

  const warned = SC.computeStepNavStatuses({
    id: 'ds', steps: ['qc', 'filter', 'pca', 'neighbors', 'cluster', 'embed'],
    warnings: ['聚类数 35 偏多（参考 20–30）'],
  })
  assert.equal(warned.cluster.warn, true)
  const qcWarn = SC.computeStepNavStatuses({ id: 'ds', steps: ['qc', 'filter'], warnings: ['线粒体比例偏高'] })
  assert.equal(qcWarn.read.warn, true)
  assert.match(SC.stepNavStatusClass('done', false), /sc-rail-done/)
  assert.match(SC.stepNavStatusClass('partial', true), /sc-rail-warn/)
})

test('buildStep2ClusterParams / EmbedParams：integration / harmony / batch_key / dims', () => {
  const harmony = SC.buildStep2ClusterParams(SC.normalizeStep2Settings({ integration: 'harmony', resolution: 0.8 }), 'sample')
  assert.equal(harmony.engine, 'seurat')
  assert.equal(harmony.integration, 'harmony')
  assert.equal(harmony.use_harmony, true)
  assert.equal(harmony.resolution, 0.8)
  assert.equal(harmony.dims, 30)
  assert.equal(harmony.n_pcs, 30)
  assert.equal(harmony.n_neighbors, 30)
  assert.equal(harmony.batch_key, 'sample')

  const rpca = SC.buildStep2ClusterParams(SC.normalizeStep2Settings({ integration: 'rpca', resolution: 0.8 }))
  assert.equal(rpca.integration, 'rpca')
  assert.equal(rpca.use_harmony, false)
  assert.equal(rpca.batch_key, undefined)

  const none = SC.buildStep2ClusterParams(SC.normalizeStep2Settings({ integration: 'none', resolution: 0.5 }))
  assert.equal(none.integration, 'none')
  assert.equal(none.resolution, 0.5)

  const builtin = SC.buildStep2ClusterParams(SC.normalizeStep2Settings({ engine: 'builtin' }))
  assert.deepEqual(builtin, { resolution: 0.8, method: 'louvain', n_top_genes: 3000 })

  const embed = SC.buildStep2EmbedParams(SC.normalizeStep2Settings({ integration: 'rpca' }), 'sample')
  assert.equal(embed.layout, 'umap')
  assert.equal(embed.integration, 'rpca')
  assert.equal(embed.engine, 'seurat')
  assert.equal(embed.batch_key, 'sample')
})

test('STEP2_PRESETS / buildStep2RecipeRequest：三套预设与 recipe id', () => {
  assert.deepEqual(SC.STEP2_PRESETS.map((preset) => preset.id), ['seurat6', 'pro-rpca', 'pro-unintegrated'])
  assert.equal(SC.buildStep2RecipeRequest('seurat6').id, 'seurat6')
  assert.equal(SC.buildStep2RecipeRequest('pro-rpca').id, 'pro-rpca')
  assert.equal(SC.buildStep2RecipeRequest('pro-unintegrated').id, 'pro-unintegrated')
  assert.equal(SC.step2PresetById('pro-rpca').integration, 'rpca')
  assert.equal(SC.step2PresetById('unknown').id, 'seurat6')
  assert.equal(SC.step2PresetById('pro-unintegrated').resolution, 0.5)
})

test('buildTrialOptions / pickTrialOption / buildTrialColorBy：三套结果选择', () => {
  const obsColumns = [
    { key: 'clusters', nUnique: 8 },
    { key: 'clusters_harmony', nUnique: 12 },
    { key: 'clusters_rpca', nUnique: 10 },
    { key: 'cell_type', nUnique: 8 },
  ]
  const options = SC.buildTrialOptions(obsColumns)
  assert.deepEqual(options.map((item) => item.id), ['harmony', 'rpca', 'unintegrated'])
  assert.equal(options[0].available, true)
  assert.equal(options[0].nClusters, 12)
  assert.equal(options[2].column, 'clusters')
  assert.equal(options[2].nClusters, 8)

  const partial = SC.buildTrialOptions([{ key: 'clusters', nUnique: 5 }])
  assert.equal(partial[0].available, false)
  assert.equal(partial[2].available, true)
  assert.equal(SC.pickTrialOption(partial, 'harmony').id, 'unintegrated')
  assert.equal(SC.pickTrialOption(SC.buildTrialOptions([]), 'harmony'), null)

  assert.deepEqual(SC.buildTrialColorBy(SC.pickTrialOption(options, 'rpca')), { kind: 'obs', key: 'clusters_rpca' })
  const rows = SC.buildTrialCompareRows(options, 'clusters_harmony')
  assert.equal(rows.length, 3)
  assert.equal(rows[0].current, true)
  assert.equal(rows[1].current, false)
})

test('buildQcCardRequest：Step1 QC / 过滤参数载荷', () => {
  const qc = SC.buildQcCardRequest('qc', { min_genes: 300, max_genes: 5000, min_counts: 500, max_pct_mt: 15, min_cells: 3 })
  assert.equal(qc.step, 'qc')
  assert.equal(qc.params.min_genes, 300)
  assert.equal(qc.params.max_pct_mt, 15)
  assert.equal(qc.params.max_genes, 5000)
  assert.equal(qc.params.min_counts, 500)
  assert.equal(qc.params.min_cells, undefined)

  const filter = SC.buildQcCardRequest('filter', { min_genes: 200, max_genes: null, min_counts: null, max_pct_mt: 20, min_cells: 3 })
  assert.equal(filter.step, 'filter')
  assert.equal(filter.params.min_genes, 200)
  assert.equal(filter.params.max_genes, undefined)
  assert.equal(filter.params.min_counts, undefined)
  assert.equal(filter.params.min_cells, 3)
  const defaults = SC.normalizeQcCardValues({})
  assert.equal(defaults.min_genes, 200)
  assert.equal(defaults.max_pct_mt, 20)
  assert.equal(defaults.max_genes, null)
})

test('formatLogEntries：运行日志格式化（时间/动作/耗时/ok）', () => {
  const entries = SC.formatLogEntries([
    { at: '2026-11-10T05:20:30.000Z', actor: 'user', action: 'step:cluster', detail: 'resolution=0.8', ok: true, ms: 1234 },
    { at: '2026-11-10T05:20:00.000Z', actor: 'agent', action: 'run:qc', ok: false, ms: 88 },
  ], 10)
  assert.equal(entries.length, 2)
  assert.ok(entries[0].time && entries[0].time !== '—')
  assert.equal(entries[0].action, 'step:cluster')
  assert.equal(entries[0].ok, true)
  assert.ok(entries[0].duration.length > 0)
  assert.equal(entries[1].ok, false)
  assert.equal(SC.formatLogEntries(null).length, 0)
  const single = SC.formatLogEntry({ action: 'x' })
  assert.equal(single.action, 'x')
  assert.equal(single.time, '—')
  assert.equal(single.ok, true)
})

test('appendGalleryFigure：画廊顺序（新图在前、同 URL 去重、上限 60）', () => {
  const first = SC.appendGalleryFigure([], { url: '/fig/a.png', title: 'A' }, { kind: 'qc' })
  assert.equal(first.length, 1)
  assert.equal(first[0].url, '/fig/a.png')
  const second = SC.appendGalleryFigure(first, { url: '/fig/b.png', title: 'B' }, { kind: 'dotplot' })
  assert.deepEqual(second.map((item) => item.url), ['/fig/b.png', '/fig/a.png'])
  const dedup = SC.appendGalleryFigure(second, { url: '/fig/a.png', title: 'A2' }, { kind: 'qc_v2' })
  assert.deepEqual(dedup.map((item) => item.url), ['/fig/a.png', '/fig/b.png'])
  assert.equal(dedup[0].title, 'A2')
  assert.equal(SC.appendGalleryFigure(first, null).length, 1)
  let many = []
  for (let index = 0; index < 70; index += 1) {
    many = SC.appendGalleryFigure(many, { url: '/fig/' + index + '.png' })
  }
  assert.equal(many.length, 60)
  assert.equal(many[0].url, '/fig/69.png')
  assert.equal(SC.formatGalleryCaption(dedup[0]), 'A2 · qc_v2')
})

test('最近数据集：upsert / 排序 / localStorage 往返', () => {
  const now = 1700000000000
  const list = SC.upsertRecentDataset([], { id: 'ds1', name: 'A', nObs: 100, nVars: 20, path: '/tmp/a' }, now)
  assert.equal(list[0].id, 'ds1')
  assert.equal(list[0].path, '/tmp/a')
  const list2 = SC.upsertRecentDataset(list, { id: 'ds2', name: 'B' }, now + 10)
  assert.deepEqual(list2.map((item) => item.id), ['ds2', 'ds1'])
  const list3 = SC.upsertRecentDataset(list2, { id: 'ds1', name: 'A2' }, now + 20)
  assert.deepEqual(list3.map((item) => item.id), ['ds1', 'ds2'])
  assert.equal(list3[0].name, 'A2')
  const storage = {
    data: {},
    getItem(key) { return this.data[key] || null },
    setItem(key, value) { this.data[key] = value },
  }
  assert.equal(SC.writeRecentDatasets(storage, list3), true)
  assert.deepEqual(SC.readRecentDatasets(storage).map((item) => item.id), ['ds1', 'ds2'])
  storage.data[SC.RECENT_STORAGE_KEY] = '{broken json'
  assert.deepEqual(SC.readRecentDatasets(storage), [])
  assert.equal(SC.normalizeRecentDatasets([{ id: 'x' }, { id: 'x' }]).length, 1)
})

test('静态：新 Shell（欢迎页/顶栏/步骤导航/日志/画廊/Step1/Step2）接线', () => {
  assert.match(CLIENT_SOURCE, /function ShellTopBar\(props\)/)
  assert.match(CLIENT_SOURCE, /function StepRail\(props\)/)
  assert.match(CLIENT_SOURCE, /function WelcomePage\(props\)/)
  assert.match(CLIENT_SOURCE, /function RunLogPanel\(props\)/)
  assert.match(CLIENT_SOURCE, /function GalleryPanel\(props\)/)
  assert.match(CLIENT_SOURCE, /function Step1Page\(props\)/)
  assert.match(CLIENT_SOURCE, /function Step2Page\(props\)/)
  assert.match(CLIENT_SOURCE, /function Step3Page\(props\)/)
  assert.match(CLIENT_SOURCE, /function Step4Page\(props\)/)
  assert.match(CLIENT_SOURCE, /h\(ShellTopBar, \{/)
  assert.match(CLIENT_SOURCE, /h\(StepRail, \{/)
  assert.match(CLIENT_SOURCE, /h\(WelcomePage, \{/)
  assert.match(CLIENT_SOURCE, /h\(RunLogPanel, \{/)
  assert.match(CLIENT_SOURCE, /h\(GalleryPanel, \{/)
  assert.match(CLIENT_SOURCE, /currentStep === 'cluster'/)
  assert.match(CLIENT_SOURCE, /setCurrentStep\(defaultStepForDataset\(activeDataset\)\)/)
  assert.match(CLIENT_SOURCE, /buildStep2RecipeRequest\(dimPreset\)/)
  assert.match(CLIENT_SOURCE, /buildTrialOptions\(obsColumns\)/)
  assert.match(CLIENT_SOURCE, /appendGalleryFigure\(previous, figure/)
  assert.match(CLIENT_SOURCE, /sc-shell-top/)
})

// ---------------------------------------------------------------------------
// 8.10 v7 Phase B：小亚群 / 注释表 / 画廊分组 / 移动端 / 进阶区
// ---------------------------------------------------------------------------

test('validateSubsetSelection / buildSubsetPayload：小亚群请求体', () => {
  assert.equal(SC.validateSubsetSelection([]).ok, false)
  assert.equal(SC.validateSubsetSelection(['', '  ']).ok, false)
  const cleaned = SC.validateSubsetSelection(['1', '1', ' 2 ', '', '3'])
  assert.equal(cleaned.ok, true)
  assert.deepEqual(cleaned.clusters, ['1', '2', '3'])
  const body = SC.buildSubsetPayload('ds1', ['0', '2', '0'], 'T-cells')
  assert.deepEqual(body, { dataset: 'ds1', clusters: ['0', '2'], name: 'T-cells' })
  const noName = SC.buildSubsetPayload('ds1', ['1'], '  ')
  assert.equal(noName.name, undefined)
  assert.equal(noName.clusters.length, 1)
})

test('normalizeSubsetResult：/dataset/subset 响应解析', () => {
  const result = SC.normalizeSubsetResult({
    ok: true,
    dataset: { id: 'ds_sub', name: 'sub-cluster', n_obs: 1478, n_vars: 800, steps: ['qc', 'filter', 'normalize', 'hvg'] },
    removed_cells: 3600,
    parent: 'ds1',
  })
  assert.equal(result.ok, true)
  assert.equal(result.datasetId, 'ds_sub')
  assert.equal(result.name, 'sub-cluster')
  assert.equal(result.removedCells, 3600)
  assert.equal(result.parent, 'ds1')
  assert.equal(result.nObs, 1478)
  assert.deepEqual(result.steps, ['qc', 'filter', 'normalize', 'hvg'])
  const live = SC.normalizeSubsetResult({ ok: true, dataset: { id: 'sub', n_obs: 1478, steps: ['qc'] }, subset_cells: 1478, kept_clusters: ['1', '2'], parent: 'ds1' })
  assert.equal(live.subsetCells, 1478)
  assert.deepEqual(live.keptClusters, ['1', '2'])
  assert.equal(live.removedCells, 0)
  const plain = SC.normalizeSubsetResult({ ok: false })
  assert.equal(plain.ok, false)
  assert.equal(plain.datasetId, '')
  assert.equal(plain.removedCells, 0)
})

test('annotationCoverage：注释表覆盖/待注释', () => {
  const rows = [
    { id: '0', label: '0', cellType: 'T cell' },
    { id: '1', label: '1', cellType: 'B cell' },
    { id: '2', label: '2', cellType: '' },
    { id: '3', label: '3', cellType: '  ' },
  ]
  const coverage = SC.annotationCoverage(rows)
  assert.equal(coverage.total, 4)
  assert.equal(coverage.mapped, 2)
  assert.deepEqual(coverage.pending, ['2', '3'])
  assert.equal(coverage.percent, 50)
  assert.equal(coverage.complete, false)
  const empty = SC.annotationCoverage([])
  assert.equal(empty.total, 0)
  assert.equal(empty.percent, 0)
  assert.equal(empty.complete, false)
  assert.equal(SC.annotationCoverage([{ id: '0', cellType: 'X' }]).complete, true)
  const zhDict = {
    cover: '已注释 {mapped}/{total} 簇（{percent}%）',
  }
  const translate = (key) => zhDict[key] || key
  assert.equal(SC.tfv(translate, 'cover', { mapped: 2, total: 4, percent: 50 }), '已注释 2/4 簇（50%）')
  assert.equal(SC.tfv(translate, 'cover', { mapped: 0 }), '已注释 0/ 簇（%）')
})

test('groupGalleryFigures：画廊按来源分组且保持首次出现顺序', () => {
  const figures = [
    { url: '/f/a.png', source: 'Step 2' },
    { url: '/f/b.png', source: 'Step 3' },
    { url: '/f/c.png', source: 'Step 2' },
    { url: '/f/d.png', source: '' },
    { url: '/f/dup.png', source: 'Step 3' },
    { url: '/f/dup.png', source: 'Step 3' },
  ]
  const groups = SC.groupGalleryFigures(figures)
  assert.deepEqual(groups.map((group) => group.source), ['Step 2', 'Step 3', ''])
  assert.deepEqual(groups[0].figures.map((figure) => figure.url), ['/f/a.png', '/f/c.png'])
  assert.equal(groups[1].figures.length, 3)
  assert.deepEqual(SC.groupGalleryFigures([]), [])
})

test('layoutModeForWidth / buildExportOptions：移动端判定与 Step4 导出顺序', () => {
  assert.equal(SC.layoutModeForWidth(1280), 'wide')
  assert.equal(SC.layoutModeForWidth(900), 'wide')
  assert.equal(SC.layoutModeForWidth(899), 'narrow')
  assert.equal(SC.layoutModeForWidth(375), 'narrow')
  assert.equal(SC.layoutModeForWidth(undefined), 'wide')
  assert.equal(SC.layoutModeForWidth('abc'), 'wide')
  assert.deepEqual(SC.EXPORT_QUICK_KINDS, ['h5ad', 'protocol', 'methods', 'report', 'script'])
  const options = SC.buildExportOptions()
  assert.deepEqual(options.map((option) => option.kind), ['h5ad', 'protocol', 'methods', 'report', 'script'])
  assert.deepEqual(options.map((option) => option.labelKey), ['exportH5ad', 'exportProtocol', 'exportMethods', 'exportReport', 'exportScript'])
})

test('静态：Phase B（小亚群 API / Step3 内联 UMAP / 进阶 PipelinePanel / 移动端）', () => {
  assert.match(CLIENT_SOURCE, /subsetDataset: function \(datasetId, clusters, name, extra\)/)
  assert.match(CLIENT_SOURCE, /jsonPost\('\/dataset\/subset'/)
  assert.match(CLIENT_SOURCE, /validateSubsetSelection\(subsetClusters\)/)
  assert.match(CLIENT_SOURCE, /setSubsetMessage\(tfv\(T, 'step3SubsetSuccess'/)
  assert.match(CLIENT_SOURCE, /handleSelect\(result\.datasetId\)/)
  assert.match(CLIENT_SOURCE, /setCurrentStep\('cluster'\)/)
  // Step3 内联 UMAP + 注释表 + 小亚群入口
  assert.match(CLIENT_SOURCE, /props\.embeddingNode \|\| h\('div', \{ className: 'sc-notice-card' \}/)
  assert.match(CLIENT_SOURCE, /T\('step3AnnotationTable'\)/)
  assert.match(CLIENT_SOURCE, /T\('step3SubsetTitle'\)/)
  // 进阶 PipelinePanel 折叠区
  assert.match(CLIENT_SOURCE, /h\(PipelinePanel, \{/)
  assert.match(CLIENT_SOURCE, /T\('advancedSteps'\)/)
  assert.match(CLIENT_SOURCE, /advancedNode: pipelineNode/)
  // 移动端单列 + 布局模式
  assert.match(CLIENT_SOURCE, /@media \(max-width:900px\)/)
  assert.match(CLIENT_SOURCE, /layoutModeForWidth\(window\.innerWidth\)/)
  assert.match(CLIENT_SOURCE, /sc-app-narrow/)
  // 旧 7-tab 死代码已清理
  assert.doesNotMatch(CLIENT_SOURCE, /function TopBar\(props\)/)
  assert.doesNotMatch(CLIENT_SOURCE, /function TabsBar\(props\)/)
  assert.doesNotMatch(CLIENT_SOURCE, /function ClustersTab\(props\)/)
  assert.doesNotMatch(CLIENT_SOURCE, /function LogTab\(props\)/)
})

// ---------------------------------------------------------------------------
// 8.11 v8：Step3 三段式（marker 方案 / 看表达 / cluster 分配）+ HVG 降权
// ---------------------------------------------------------------------------

test('normalizeSchemeGenes / schemeRowsFromPreset：基因解析与内置方案', () => {
  assert.deepEqual(SC.normalizeSchemeGenes('CD3D, CD3E; IL7R\nNKG7   GNLY'), ['CD3D', 'CD3E', 'IL7R', 'NKG7', 'GNLY'])
  assert.deepEqual(SC.normalizeSchemeGenes('cd3d CD3D Cd3d'), ['cd3d'])
  assert.deepEqual(SC.normalizeSchemeGenes('CD3D cd3d Cd3D'), ['CD3D'], '大小写不敏感去重且保留首个写法')
  assert.deepEqual(SC.normalizeSchemeGenes(''), [])
  assert.deepEqual(SC.normalizeSchemeGenes(['A', 'a', 'B']), ['A', 'B'])

  const human = SC.MARKER_PRESETS.human.rows
  const mouse = SC.MARKER_PRESETS.mouse.rows
  assert.ok(human.length >= 18, '人方案 ≥ 18 类')
  assert.ok(mouse.length >= 18, '鼠方案 ≥ 18 类')
  for (const name of ['Mast cell', 'Neutrophil', 'DC', 'Monocyte', 'Macrophage', 'B cell', 'Plasma cell', 'T cell', 'NK cell', 'Epithelial', 'Endothelial cells', 'CAF', 'Smooth muscle cells', 'Tumor/Epithelial']) {
    assert.ok(human.some((row) => row.cellType === name), '人方案缺少：' + name)
  }
  assert.ok(mouse.some((row) => row.cellType === 'Mast cell'))
  assert.ok(mouse.some((row) => row.cellType === 'Tumor cell'))
  const humanRows = SC.schemeRowsFromPreset('human')
  assert.equal(humanRows.length, human.length)
  assert.ok(humanRows[0].genesText.length > 0)
  assert.equal(SC.schemeRowsFromPreset('unknown').length, human.length)
})

test('buildMarkerScheme：校验（空类型/重复类型/无基因/超 50 行）', () => {
  const scheme = SC.buildMarkerScheme([
    { cellType: 'T cell', genesText: 'CD3D CD3E' },
    { cellType: 'B cell', genesText: 'MS4A1 CD79A' },
  ])
  assert.equal(scheme.valid, true)
  assert.equal(scheme.errors.length, 0)
  assert.equal(scheme.geneCount, 4)
  assert.equal(scheme.typedRows, 2)
  assert.equal(scheme.usableRows.length, 2)

  const bad = SC.buildMarkerScheme([
    { cellType: '', genesText: 'CD3D' },
    { cellType: 'T cell', genesText: 'CD3D' },
    { cellType: 't cell', genesText: 'CD3E' },
    { cellType: 'B cell', genesText: '' },
  ])
  assert.equal(bad.valid, false)
  assert.ok(bad.errors.some((item) => item === 'EMPTY_TYPE'))
  assert.ok(bad.errors.some((item) => item.indexOf('DUPLICATE:') === 0))
  assert.ok(bad.warnings.some((item) => item === 'NO_GENES:B cell'))

  const many = []
  for (let index = 0; index < 55; index += 1) many.push({ cellType: 'T' + index, genesText: 'CD3D' })
  assert.ok(SC.buildMarkerScheme(many).warnings.includes('TOO_MANY_ROWS'))
  assert.equal(SC.buildMarkerScheme([]).geneCount, 0)
  assert.equal(SC.SCHEME_MAX_ROWS, 50)
})

test('flattenSchemeGenes / buildGroupedDotplotParams：扁平化与分组出图参数', () => {
  const scheme = SC.buildMarkerScheme([
    { cellType: 'T cell', genesText: 'CD3D CD3E' },
    { cellType: 'B cell', genesText: 'cd3d MS4A1' },
  ])
  assert.deepEqual(SC.flattenSchemeGenes(scheme), ['CD3D', 'CD3E', 'MS4A1'], '跨类型全局去重')
  const params = SC.buildGroupedDotplotParams(scheme)
  assert.equal(params.groupby, 'clusters')
  assert.deepEqual(params.genes, ['CD3D', 'CD3E', 'MS4A1'])
  assert.deepEqual(params.gene_groups, [
    { label: 'T cell', genes: ['CD3D', 'CD3E'] },
    { label: 'B cell', genes: ['cd3d', 'MS4A1'] },
  ])
  const empty = SC.buildGroupedDotplotParams(SC.buildMarkerScheme([{ cellType: 'X', genesText: '' }]))
  assert.deepEqual(empty.genes, [])
  assert.deepEqual(empty.gene_groups, [])
})

test('normalizeMarkerCsv：cell_type,genes 两列解析', () => {
  const parsed = SC.normalizeMarkerCsv('cell_type,genes\nT cell,CD3D CD3E\nB cell,MS4A1\n')
  assert.equal(parsed.rows.length, 2)
  assert.deepEqual(parsed.rows[0], { cellType: 'T cell', genesText: 'CD3D CD3E' })
  const noHeader = SC.normalizeMarkerCsv('T cell,CD3D\nB cell,MS4A1 CD79A')
  assert.equal(noHeader.rows.length, 2)
  const tabs = SC.normalizeMarkerCsv('NK cell\tNKG7,GNLY\tPRF1')
  assert.equal(tabs.rows[0].cellType, 'NK cell')
  assert.equal(tabs.rows[0].genesText, 'NKG7 GNLY PRF1')
  assert.deepEqual(SC.normalizeMarkerCsv('').rows, [])
  assert.deepEqual(SC.normalizeMarkerCsv(',CD3D').errors, ['EMPTY_TYPE'])
})

test('buildClusterAssignmentRows / buildClusterAnnotation：cluster → 类型映射', () => {
  const scheme = SC.buildMarkerScheme([
    { cellType: 'T cell', genesText: 'CD3D' },
    { cellType: 'B cell', genesText: 'MS4A1' },
  ])
  const clusterRows = [
    { id: '0', label: '0', nCells: 312, topGenes: ['CD3D', 'CD3E'] },
    { id: '1', label: '1', nCells: 120, topGenes: [] },
  ]
  const rows = SC.buildClusterAssignmentRows(clusterRows, scheme, { 0: 'T cell', 1: '__custom__' }, { 1: 'Myoepithelial' })
  assert.deepEqual(rows[0].options, ['T cell', 'B cell'])
  assert.equal(rows[0].assignment, 'T cell')
  assert.equal(rows[1].assignment, '__custom__')
  assert.equal(rows[1].custom, 'Myoepithelial')
  assert.deepEqual(rows[0].topGenes, ['CD3D', 'CD3E'])

  assert.deepEqual(SC.buildClusterAnnotation({ 0: 'T cell', 1: '__custom__' }, { 1: 'Myoepithelial' }), { 0: 'T cell', 1: 'Myoepithelial' })
  assert.deepEqual(SC.buildClusterAnnotation({ 0: '', 1: '__custom__', 2: 'B cell' }, {}), { 2: 'B cell' })
  assert.deepEqual(SC.buildClusterAnnotation({}, {}), {})
})

test('HVG 参考参数：defaultDimredParams / buildDimredPayload 仍带 3000', () => {
  assert.equal(SC.defaultDimredParams().nTopGenes, 3000)
  assert.equal(SC.STEP2_DEFAULTS.nTopGenes, 3000)
  const payload = SC.buildDimredPayload(SC.defaultDimredParams(), 'sample')
  assert.equal(payload.n_top_genes, 3000)
  assert.equal(payload.engine, 'seurat')
  assert.equal(payload.integration, 'harmony')
  assert.equal(payload.use_harmony, true)
  assert.equal(payload.resolution, 0.8)
  assert.equal(payload.n_pcs, 30)
  assert.equal(payload.dims, 30)
  assert.equal(payload.n_neighbors, 30)
  assert.equal(payload.batch_key, 'sample')
  // 折叠 UI 改值后仍会生效
  const custom = SC.buildDimredPayload(SC.normalizeStep2Settings({ nTopGenes: 2000, engine: 'builtin' }))
  assert.equal(custom.n_top_genes, 2000)
  assert.equal(custom.engine, 'builtin')
  assert.equal(custom.use_harmony, undefined)
  // 只重聚类 payload 也带 n_top_genes
  assert.equal(SC.buildStep2ClusterParams(SC.defaultDimredParams()).n_top_genes, 3000)
  assert.equal(SC.buildStep2ClusterParams(SC.normalizeStep2Settings({ engine: 'builtin' })).n_top_genes, 3000)
})

test('静态：Step2 HVG 降权（主参数区无 n_top_genes，折叠区保留 hint）', () => {
  const start = CLIENT_SOURCE.indexOf('function Step2Page(props) {')
  const end = CLIENT_SOURCE.indexOf('function Step3Page(props) {', start)
  assert.ok(start > 0 && end > start, 'Step2Page 代码段存在')
  const snippet = CLIENT_SOURCE.slice(start, end)
  const advancedIndex = snippet.indexOf("T('step2AdvancedParams')")
  const topGenesIndex = snippet.indexOf("T('dimredTopGenes')")
  assert.ok(advancedIndex > 0, '进阶折叠区存在')
  assert.ok(topGenesIndex > advancedIndex, 'n_top_genes 只在折叠区出现（排在进阶标题之后）')
  assert.ok(snippet.indexOf("T('step2HvgHint')") > advancedIndex, '折叠区有 hint 文案')
  assert.match(SC.DICT_ZH.step2HvgHint, /3000/)
  assert.ok(snippet.indexOf('hvgReference') > 0, '页面有 hvg 参考行')
  // 主区字段只剩 引擎/整合/预设/batch/resolution/n_pcs/neighbors
  const mainTopGenes = snippet.slice(0, advancedIndex).indexOf("T('dimredTopGenes')")
  assert.equal(mainTopGenes, -1, '主参数区不应出现 n_top_genes')
})

test('静态：Step3 三段式顺序（方案卡 → 看表达 → cluster 分配）与 gene_groups 接线', () => {
  const start = CLIENT_SOURCE.indexOf('function Step3Page(props) {')
  const end = CLIENT_SOURCE.indexOf('function Step4Page(props) {', start)
  const snippet = CLIENT_SOURCE.slice(start, end)
  const schemeIndex = snippet.indexOf('h(MarkerSchemeCard, {')
  const expressionIndex = snippet.indexOf('h(ExpressionCard, {')
  const assignIndex = snippet.indexOf('h(ClusterAssignCard, {')
  assert.ok(schemeIndex >= 0 && expressionIndex > schemeIndex && assignIndex > expressionIndex, '三段顺序：方案 → 看表达 → 分配')
  assert.ok(snippet.indexOf('h(AnnotationWorkbench') < 0, '旧「先列 cluster」工作台不再出现在 Step3 主流程')
  assert.match(CLIENT_SOURCE, /buildGroupedDotplotParams\(built\)/)
  assert.match(CLIENT_SOURCE, /gene_groups/)
  assert.match(CLIENT_SOURCE, /handleSchemeChange/)
  assert.match(CLIENT_SOURCE, /handleSaveClusterAssignments/)
  assert.match(CLIENT_SOURCE, /buildMarkerScheme\(schemeRows\)/)
  assert.match(CLIENT_SOURCE, /MARKER_PRESETS/)
})

// ---------------------------------------------------------------------------
// 9. 静态合规：模块形态 / 槽注册 / 样式约束
// ---------------------------------------------------------------------------

test('静态：ModuleLoader 形态、只 require react、无 JSX、无 @deepseek-ai 依赖', () => {
  assert.match(CLIENT_SOURCE, /window\.__ModuleLoader__\.load\(\{ id: NS, factory: factory \}\)/)
  assert.match(CLIENT_SOURCE, /var NS = 'dsh-singlecell-workbench'/)
  assert.match(CLIENT_SOURCE, /var API = '\/singlecell'/)
  assert.match(CLIENT_SOURCE, /factory\(require\)/)
  assert.match(CLIENT_SOURCE, /var React = require\('react'\)/)
  assert.doesNotMatch(CLIENT_SOURCE, /@deepseek-ai\//)
  // 只允许 require('react')：检查所有 require 调用
  const requires = [...CLIENT_SOURCE.matchAll(/require\('([^']+)'\)/g)].map((match) => match[1])
  assert.deepEqual([...new Set(requires)], ['react'])
  // 无 JSX：不出现 <div>/<span> 这种在 JS 里裸写的标签
  assert.doesNotMatch(CLIENT_SOURCE, /return\s*<\s*[A-Za-z]/)
})

test('静态：槽注册 main(key=singlecell) + sidebar.panellist(id=singlecell, order=3) + locale/effect 清理', () => {
  assert.match(CLIENT_SOURCE, /inject: \['slots', 'locale'\]/)
  assert.match(CLIENT_SOURCE, /slots\.inject\('main'/)
  assert.match(CLIENT_SOURCE, /name: 'main', key: 'singlecell', locale: NS/)
  assert.match(CLIENT_SOURCE, /slots\.inject\('sidebar\.panellist'/)
  assert.match(CLIENT_SOURCE, /name: 'sidebar\.panellist',\s*\n\s*id: 'singlecell',\s*\n\s*order: 3,/)
  assert.match(CLIENT_SOURCE, /locale\.register\(NS, \{ zh: DICT_ZH, en: DICT_EN \}\)/)
  assert.match(CLIENT_SOURCE, /ctx\.effect\(registerDictionaries, 'dsh-singlecell-workbench: locale'\)/)
  assert.match(CLIENT_SOURCE, /ctx\.effect\(function \(\) \{/)
  assert.match(CLIENT_SOURCE, /PluginErrorBoundary/)
  assert.match(CLIENT_SOURCE, /getDerivedStateFromError/)
})

test('静态：三栏布局 + 主题变量 + sc- 前缀 + 无硬编码颜色', () => {
  assert.match(CLIENT_SOURCE, /height:100%;min-height:72vh/)
  assert.match(CLIENT_SOURCE, /--dsw-alias-brand-primary/)
  assert.match(CLIENT_SOURCE, /--dsw-alias-bg-layer-1/)
  assert.match(CLIENT_SOURCE, /--dsw-alias-border-l1/)
  assert.match(CLIENT_SOURCE, /--dsw-alias-state-error-primary/)
  assert.match(CLIENT_SOURCE, /--dsw-alias-state-warn-primary/)
  assert.match(CLIENT_SOURCE, /color-mix\(in srgb,var\(--dsw-alias-brand-primary\)/)
  assert.match(CLIENT_SOURCE, /flex:0 0 268px/)
  assert.match(CLIENT_SOURCE, /flex:0 0 300px/)

  const cssStart = CLIENT_SOURCE.indexOf('var CSS = [')
  const cssEnd = CLIENT_SOURCE.indexOf("].join('')", cssStart)
  assert.ok(cssStart > 0 && cssEnd > cssStart, '应能定位 CSS 字符串')
  const cssBlock = CLIENT_SOURCE.slice(cssStart, cssEnd)
  const selectors = [...cssBlock.matchAll(/'\.([A-Za-z][\w-]*)/g)].map((match) => match[1])
  assert.ok(selectors.length > 60, `CSS 选择器太少：${selectors.length}`)
  for (const selector of selectors) {
    assert.ok(selector.startsWith('sc-'), `类名未用 sc- 前缀：${selector}`)
  }
  assert.doesNotMatch(cssBlock, /#[0-9a-fA-F]{3,6}/, 'CSS 里不应出现硬编码 #hex 颜色')
})

test('静态：交互与轮询要求存在（canvas / visibility / interval 清理 / clipboard / 7 tab）', () => {
  assert.match(CLIENT_SOURCE, /getContext\('2d'\)/)
  assert.match(CLIENT_SOURCE, /devicePixelRatio/)
  assert.match(CLIENT_SOURCE, /addEventListener\('wheel', handleWheel, \{ passive: false \}\)/)
  assert.match(CLIENT_SOURCE, /onDoubleClick/)
  assert.match(CLIENT_SOURCE, /document\.visibilityState === 'hidden'/)
  assert.match(CLIENT_SOURCE, /setInterval\(tick, POLL_MS\)/)
  assert.match(CLIENT_SOURCE, /clearInterval\(timer\)/)
  assert.match(CLIENT_SOURCE, /visibilitychange/)
  assert.match(CLIENT_SOURCE, /navigator\.clipboard/)
  assert.match(CLIENT_SOURCE, /x-dsh-filename/)
  assert.match(CLIENT_SOURCE, /accept: IMPORT_ACCEPT/)
  // 7 个 tab（v2 六个 + v3 高级分析）
  for (const key of ['tabEmbedding', 'tabQc', 'tabMarkers', 'tabClusters', 'tabRecord', 'tabLog', 'tabAnalysis']) {
    assert.ok(CLIENT_SOURCE.includes("labelKey: '" + key + "'"), `缺少 tab：${key}`)
  }
  // v2 + v3 路由都在客户端接到
  for (const path of ['/audit', '/autoannotate', '/resolution-scan', '/recipe', '/formats', '/import/inspect', '/analyses', '/analysis']) {
    assert.ok(CLIENT_SOURCE.includes("'" + path + "'"), `缺少路由：${path}`)
  }
  assert.doesNotMatch(CLIENT_SOURCE, /document\.body\.appendChild/)
  assert.doesNotMatch(CLIENT_SOURCE, /\balert\(/)
})

test('静态：v3 导入弹层 + 高级分析 tab 的接线', () => {
  assert.match(CLIENT_SOURCE, /function ImportDialog\(props\)/)
  assert.match(CLIENT_SOURCE, /function AnalysisTab\(props\)/)
  assert.match(CLIENT_SOURCE, /className: 'sc-modal'/)
  assert.match(CLIENT_SOURCE, /role: 'dialog'/)
  assert.match(CLIENT_SOURCE, /'aria-modal': 'true'/)
  assert.match(CLIENT_SOURCE, /function handleImport\(options\)/)
  assert.match(CLIENT_SOURCE, /function handlePickCandidate\(candidate\)/)
  assert.match(CLIENT_SOURCE, /function handleRunAnalysis\(\)/)
  assert.match(CLIENT_SOURCE, /colorOptions: colorOptions/)
  // 只 require('react')（不使用 react-dom portal），弹层用 fixed 浮层
  assert.doesNotMatch(CLIENT_SOURCE, /require\('react-dom'\)/)
  assert.match(CLIENT_SOURCE, /\.sc-modal\{position:fixed;inset:0;/)
})
