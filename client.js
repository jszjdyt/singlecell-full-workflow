/**
 * dsh-singlecell — 浏览器端 UI（Client 半边）
 *
 * 形态：window.__ModuleLoader__.load({ id:'dsh-singlecell', factory(require){...} })，
 * 只 require('react')，不碰任何 @deepseek-ai Client 包（宿主版本一变就会炸）。
 *
 * 数据面：全部走同源相对路径 /singlecell/*（见 CONTRACT.md §3 / §6）。
 *
 * 文件分两层：
 *   1) 纯函数层（本段）：STEPS 定义 / 数值与时间格式化 / 坐标变换 / 归一化 /
 *      建议→中文指令 / 可注入 fetch 的 API 客户端……不碰 DOM、不 require react。
 *      node 环境（typeof window === 'undefined'）挂到 globalThis.__SC_TEST__，
 *      供 test/client-logic.test.mjs 直接 import 测试。
 *   2) 浏览器层（factory 内）：React 组件 + 样式 + 槽注册。顶层的
 *      `if (typeof window === 'undefined')` 守卫保证 node 能安全 import 本文件。
 *
 * 约定：所有类名 sc- 前缀；颜色只用 --dsw-alias-* 主题变量（见 CSS 段），
 * 服务端下发的聚类色（categories[].color）与 JS 常量 PALETTE 属于「数据色」例外。
 */
;(function () {
  'use strict'

  var NS = 'dsh-singlecell'
  var API = '/singlecell'
  var POLL_MS = 2500
  var DEFAULT_API_TIMEOUT_MS = 60000
  var HEAVY_API_TIMEOUT_MS = 300000
  var MAX_HISTORY_ROWS = 20

  /** 聚类调色板（数据色）：服务端会下发 categories[].color，这里只是缺省兜底。 */
  var PALETTE = [
    '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948',
    '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac', '#1f77b4', '#ff7f0e',
    '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf',
  ]

  /** 连续着色色带（数据色）：低→高。 */
  var CONTINUOUS_STOPS = ['#3b4cc0', '#6a9fdb', '#9ec6e8', '#e8d4a8', '#f0a35e', '#c0392b']

  /** 自检结果的三种等级样式（warn 黄 / fail 红，由主题色派生）。 */
  var AUDIT_LEVEL_STYLE = {
    ok: { tone: 'ok', icon: '✓', labelKey: 'auditOk' },
    warn: { tone: 'warn', icon: '!', labelKey: 'auditWarn' },
    fail: { tone: 'fail', icon: '×', labelKey: 'auditFail' },
  }

  /** 记录 tab 的导出按钮组（顺序固定，与 CONTRACT §6.2 的 export what 对应）。 */
  var EXPORT_KINDS = [
    { id: 'protocol', labelKey: 'exportProtocol', hint: 'JSON' },
    { id: 'methods', labelKey: 'exportMethods', hint: 'Markdown' },
    { id: 'report', labelKey: 'exportReport', hint: 'Markdown' },
    { id: 'script', labelKey: 'exportScript', hint: 'Python' },
    { id: 'h5ad', labelKey: 'exportH5ad', hint: 'h5ad' },
    { id: 'markers_csv', labelKey: 'exportMarkersCsv', hint: 'CSV' },
  ]

  // ==========================================================================
  // 流程步骤：顺序即展示顺序；requires 即依赖（与 CONTRACT §2.2 表逐字一致）。
  // params 里 camel/snake 的 key 必须与 Python worker 的字段一致。
  // ==========================================================================
  var STEPS = [
    {
      id: 'qc',
      titleKey: 'stepQc',
      title: '质控 QC',
      requires: [],
      params: [
        { key: 'min_genes', kind: 'number', labelKey: 'paramMinGenes', default: 200, min: 0, step: 10 },
        { key: 'max_genes', kind: 'number', labelKey: 'paramMaxGenes', default: null, min: 0, step: 100 },
        { key: 'min_counts', kind: 'number', labelKey: 'paramMinCounts', default: null, min: 0, step: 100 },
        { key: 'max_pct_mt', kind: 'number', labelKey: 'paramMaxPctMt', default: 20, min: 0, max: 100, step: 1 },
      ],
    },
    {
      id: 'filter',
      titleKey: 'stepFilter',
      title: '过滤',
      requires: ['qc'],
      params: [
        { key: 'min_genes', kind: 'number', labelKey: 'paramMinGenes', default: 200, min: 0, step: 10 },
        { key: 'max_genes', kind: 'number', labelKey: 'paramMaxGenes', default: null, min: 0, step: 100 },
        { key: 'min_counts', kind: 'number', labelKey: 'paramMinCounts', default: null, min: 0, step: 100 },
        { key: 'max_pct_mt', kind: 'number', labelKey: 'paramMaxPctMt', default: 20, min: 0, max: 100, step: 1 },
        { key: 'min_cells', kind: 'number', labelKey: 'paramMinCells', default: 3, min: 0, step: 1 },
      ],
    },
    {
      id: 'normalize',
      titleKey: 'stepNormalize',
      title: '标准化',
      requires: ['filter'],
      params: [
        { key: 'target_sum', kind: 'number', labelKey: 'paramTargetSum', default: 10000, min: 1, step: 1000 },
      ],
    },
    {
      id: 'hvg',
      titleKey: 'stepHvg',
      title: '高变基因',
      requires: ['normalize'],
      params: [
        { key: 'n_top_genes', kind: 'number', labelKey: 'paramNTopGenes', default: 2000, min: 10, step: 100 },
        {
          key: 'flavor',
          kind: 'select',
          labelKey: 'paramFlavor',
          default: 'seurat',
          options: [
            { value: 'seurat', labelKey: 'flavorSeurat' },
            { value: 'cell_ranger', labelKey: 'flavorCellRanger' },
          ],
        },
      ],
    },
    {
      id: 'pca',
      titleKey: 'stepPca',
      title: 'PCA',
      requires: ['hvg'],
      params: [
        { key: 'n_pcs', kind: 'number', labelKey: 'paramNPcs', default: 30, min: 2, step: 1 },
      ],
    },
    {
      id: 'neighbors',
      titleKey: 'stepNeighbors',
      title: '邻接图',
      requires: ['pca'],
      params: [
        { key: 'n_neighbors', kind: 'number', labelKey: 'paramNNeighbors', default: 15, min: 2, step: 1 },
        {
          key: 'metric',
          kind: 'select',
          labelKey: 'paramMetric',
          default: 'cosine',
          options: [
            { value: 'cosine', labelKey: 'metricCosine' },
            { value: 'euclidean', labelKey: 'metricEuclidean' },
            { value: 'manhattan', labelKey: 'metricManhattan' },
          ],
        },
      ],
    },
    {
      id: 'cluster',
      titleKey: 'stepCluster',
      title: '聚类',
      requires: ['neighbors'],
      params: [
        { key: 'resolution', kind: 'number', labelKey: 'paramResolution', default: 1.0, min: 0.05, max: 5, step: 0.05 },
        {
          key: 'method',
          kind: 'select',
          labelKey: 'paramMethod',
          default: 'louvain',
          options: [
            { value: 'louvain', labelKey: 'methodLouvain' },
            { value: 'leiden', labelKey: 'methodLeiden' },
          ],
        },
      ],
    },
    {
      id: 'embed',
      titleKey: 'stepEmbed',
      title: '二维嵌入',
      requires: ['pca'],
      params: [
        {
          key: 'layout',
          kind: 'select',
          labelKey: 'paramLayout',
          default: 'tsne',
          options: [
            { value: 'tsne', labelKey: 'layoutTsne' },
            { value: 'pca', labelKey: 'layoutPca' },
            { value: 'umap', labelKey: 'layoutUmap' },
          ],
        },
        { key: 'perplexity', kind: 'number', labelKey: 'paramPerplexity', default: 30, min: 2, step: 1 },
      ],
    },
    {
      id: 'markers',
      titleKey: 'stepMarkers',
      title: '标记基因',
      requires: ['cluster'],
      params: [
        { key: 'groupby', kind: 'select', labelKey: 'paramGroupby', default: 'clusters', options: [{ value: 'clusters', labelKey: 'groupbyClusters' }] },
        { key: 'top_n', kind: 'number', labelKey: 'paramTopN', default: 10, min: 1, max: 50, step: 1 },
        {
          key: 'method',
          kind: 'select',
          labelKey: 'paramMarkerMethod',
          default: 't-test',
          options: [
            { value: 't-test', labelKey: 'markerTTest' },
            { value: 'wilcoxon', labelKey: 'markerWilcoxon' },
          ],
        },
      ],
    },
    {
      id: 'annotate',
      titleKey: 'stepAnnotate',
      title: '细胞类型注释',
      requires: ['cluster'],
      params: [],
    },
  ]

  var STEP_BY_ID = {}
  for (var stepIndex = 0; stepIndex < STEPS.length; stepIndex++) STEP_BY_ID[STEPS[stepIndex].id] = STEPS[stepIndex]

  /** 纯函数层用的中文步骤名（buildAgentInstruction / history / recipe 进度）。 */
  var STEP_TITLES_ZH = {
    qc: '质控',
    filter: '过滤',
    normalize: '标准化',
    hvg: '高变基因',
    pca: 'PCA',
    neighbors: '邻接图',
    cluster: '聚类',
    embed: '二维嵌入',
    markers: '标记基因',
    annotate: '细胞类型注释',
  }

  // ==========================================================================
  // 语言字典：中文为主，英文兜底。可见文案全部经 T()/tf() 取。
  // 占位符统一 `{v}`（由 tf() 替换；宿主 t() 的插值行为不可依赖）。
  // ==========================================================================
  var DICT_ZH = {
    nav: '单细胞',
    title: '单细胞分析工作台',
    subtitle: 'QC · 聚类 · 标记基因 · 嵌入图',
    emptyTitle: '还没有数据集',
    emptyHint: '载入内置示例数据立刻试一遍，或输入 .h5ad 绝对路径打开自己的数据。',
    loadDemo: '载入示例数据',
    openH5ad: '打开 h5ad',
    pathPlaceholder: '/Users/you/data/sample.h5ad',
    refresh: '刷新',
    loading: '加载中…',
    opening: '打开中…',
    uploading: '上传中…',
    openPath: '打开',
    uploadFile: '上传文件…',
    engine: '引擎',
    engineBuiltin: '内置引擎',
    engineScanpy: 'scanpy',
    engineNone: 'Python 环境未就绪',
    engineNotReady: '未就绪',
    pipeline: '分析流程',
    pipelineHint: '点卡片标题展开参数；已完成的步骤可以重跑。',
    agentPanel: 'Agent',
    agentHint: '也可以在对话里直接让 Agent 操作这个工作台。',
    suggestions: '建议下一步',
    copyInstruction: '复制指令',
    copied: '已复制',
    activity: 'Agent 活动',
    tabEmbedding: '嵌入图',
    tabQc: 'QC',
    tabMarkers: '标记基因',
    tabClusters: '聚类',
    tabRecord: '记录',
    tabLog: '日志',
    colorBy: '着色',
    run: '运行',
    rerun: '重跑',
    running: '运行中…',
    done: '已完成',
    pending: '待运行',
    needPrereq: '需前置',
    generateQc: '生成 QC 图',
    generateDotplot: '生成 dot plot',
    cellTypePlaceholder: '细胞类型',
    noData: '先载入数据',
    // 步骤
    stepQc: '质控 QC',
    stepFilter: '过滤细胞',
    stepNormalize: '标准化',
    stepHvg: '高变基因',
    stepPca: 'PCA 降维',
    stepNeighbors: '邻接图',
    stepCluster: '聚类',
    stepEmbed: '二维嵌入',
    stepMarkers: '标记基因',
    stepAnnotate: '细胞类型注释',
    // 参数
    paramMinGenes: '最少基因数',
    paramMaxGenes: '最多基因数',
    paramMinCounts: '最少 UMI',
    paramMaxPctMt: '线粒体占比上限 (%)',
    paramMinCells: '基因最少细胞数',
    paramTargetSum: '目标文库大小',
    paramNTopGenes: '高变基因数',
    paramFlavor: 'HVG 算法',
    paramNPcs: '主成分数',
    paramNNeighbors: '近邻数',
    paramMetric: '距离度量',
    paramResolution: '聚类分辨率',
    paramMethod: '聚类方法',
    paramLayout: '嵌入布局',
    paramPerplexity: '困惑度',
    paramGroupby: '分组列',
    paramTopN: '每组基因数',
    paramMarkerMethod: '检验方法',
    flavorSeurat: 'seurat',
    flavorCellRanger: 'cell_ranger',
    metricCosine: 'cosine',
    metricEuclidean: 'euclidean',
    metricManhattan: 'manhattan',
    methodLouvain: 'louvain',
    methodLeiden: 'leiden',
    layoutTsne: 't-SNE',
    layoutPca: 'PCA',
    layoutUmap: 'UMAP',
    groupbyClusters: 'clusters',
    markerTTest: 't-test',
    markerWilcoxon: 'wilcoxon',
    // 嵌入图
    noEmbeddingTitle: '还没有嵌入坐标',
    noEmbeddingHint: '运行「二维嵌入」步骤后，这里会出现可缩放、可悬停、可高亮的散点图。',
    runEmbed: '运行二维嵌入',
    embedNeedPca: '需先完成 PCA',
    embedFailed: '嵌入图加载失败',
    embedEmpty: '没有可显示的点',
    loadingEmbedding: '正在加载嵌入图…',
    downsampleNotice: '已下采样显示 {v} 个点',
    resetView: '复位',
    colorCluster: '聚类 clusters',
    colorCellType: '细胞类型 cell_type',
    colorNGenes: '基因数 n_genes',
    colorPctMt: '线粒体占比 pct_counts_mt',
    colorGene: '基因表达',
    genePlaceholder: '输入基因名，如 CD3D',
    geneEmptyHint: '输入基因名后按回车或点「着色」，按表达量给细胞着色。',
    geneApply: '着色',
    legend: '图例',
    pointsUnit: '个细胞',
    hoverHint: '滚轮缩放 · 拖拽平移 · 双击复位 · 点图例高亮',
    // QC
    qcMetrics: '指标',
    qcEmptyHint: '勾选指标后点「生成 QC 图」；图按需生成，不会自动跑。',
    qcNeedStep: '提示：先在左栏运行「质控 QC」步骤。',
    qcFigureAlt: 'QC 小提琴图',
    plotFailed: '出图失败',
    // 标记基因
    markersNeedCluster: '先运行「聚类」步骤，才能计算标记基因。',
    markersEmpty: '还没有标记基因结果',
    loadMarkers: '计算标记基因',
    reloadMarkers: '重新计算',
    dotplotAlt: '标记基因 dot plot',
    tableCluster: '簇',
    tableGene: '基因',
    tableLogfc: 'logFC',
    tablePval: 'adj.P',
    tablePctIn: 'pct_in',
    // 聚类
    clusterBarTitle: '各簇细胞数',
    generateClusterBar: '生成聚类柱状图',
    clusterBarFailed: '聚类柱状图生成失败',
    annotationTitle: '细胞类型注释',
    annotationHint: '改完细胞类型名后点保存，会写回当前数据集的 cell_type。',
    annotationNeedMapping: '请先在「聚类」页填写至少一个细胞类型名',
    annotationCountHint: '已填 {v} 个细胞类型（在「聚类」页编辑）',
    saveAnnotation: '保存注释',
    annotationSaved: '已保存',
    cellsUnit: '细胞',
    genesUnit: '基因',
    // 日志 / Agent
    logEmpty: '还没有活动记录',
    noSuggestions: '暂无建议',
    activityEmpty: '暂无活动',
    progress: '已完成',
    statusTitle: '当前状态',
    colTime: '时间',
    colStep: '步骤',
    colParams: '参数',
    colSummary: '结果',
    colMs: '耗时',
    resultOk: '成功',
    resultFail: '失败',
    copyFailed: '复制失败（浏览器不允许写剪贴板）',
    requestFailed: '请求失败',
    pollFailed: '状态刷新失败',
    errorTitle: '出错了',
    retry: '重试',
    crashTitle: '工作台显示异常',
    crashHint: '单块渲染出错，点「重试」恢复；数据集和运行记录不受影响。',
    // v2：预设流程
    recipes: '预设流程',
    recipePick: '选择预设',
    runRecipe: '一键运行',
    recipeRunning: '运行中…',
    recipeStepOf: '第 {v} 步',
    recipeDone: '预设流程已跑完 {v} 步',
    recipeFailed: '预设流程执行失败',
    // v2：记录 tab
    historyTitle: '数据集历史',
    historyEmpty: '还没有分析记录；跑一个步骤后，这里会出现可追溯的溯源历史。',
    auditTitle: '自检',
    runAudit: '运行自检',
    auditRunning: '自检中…',
    auditEmpty: '还没有自检结果',
    auditScore: '得分 {v}',
    auditOk: '通过',
    auditWarn: '注意',
    auditFail: '问题',
    autoannotateTitle: '自动注释建议',
    runAutoannotate: '生成注释建议',
    autoannotateRunning: '生成中…',
    autoannotateEmpty: '还没有注释建议',
    autoannotatePanel: '参考面板',
    autoannotateScore: 'score',
    autoannotateMargin: 'margin',
    autoannotateEvidence: '证据基因',
    autoannotateAlternative: '次选',
    adopt: '采用',
    adoptAll: '全部采用',
    adopted: '已采用',
    adoptNeedCluster: '先运行「聚类」步骤后再自动注释。',
    resolutionTitle: '分辨率扫描',
    resolutionPlaceholder: '0.4,0.8,1.2',
    runScan: '开始扫描',
    scanRunning: '扫描中…',
    scanEmpty: '还没有扫描结果',
    scanHint: '在既有 kNN 图上对比不同 resolution 的簇数与 ARI，不改变当前聚类。',
    scanCurrent: '当前',
    colResolution: 'resolution',
    colClusters: '簇数',
    colAri: 'ARI',
    exportTitle: '导出产物',
    exportProtocol: '可复现协议',
    exportMethods: 'Methods',
    exportReport: '分析报告',
    exportScript: '复现脚本',
    exportH5ad: 'h5ad 数据',
    exportMarkersCsv: '标记基因 CSV',
    exportRunning: '导出中…',
    exportDone: '已导出：',
    copyPath: '复制路径',
    recordNeedDataset: '先载入数据，记录会出现在这里。',
    // 其它
    h5adNeeded: '请选择 .h5ad / .h5 文件',
    pollErrorPrefix: '状态刷新失败',
    // v3：导入面板
    importData: '导入数据',
    importTitle: '导入单细胞数据',
    importSubtitle: 'h5ad / 10x / loom / 矩阵文本 / GEO series / R (.rds/.RData) / 压缩包 / URL',
    importSourcePath: '本地路径',
    importSourceUrl: '下载 URL',
    importSourceFile: '本地文件',
    importPathPlaceholder: '/Users/you/data/sample.h5ad 或目录',
    importPathHint: '支持文件或目录；目录会找 10x 三件套。相对路径按 Host 工作目录解析。',
    importUrlPlaceholder: 'https://ftp.ncbi.nlm.nih.gov/geo/series/.../suppl/xxx.txt.gz',
    importUrlHint: '支持 http(s)/ftp 直链（GEO 补充文件），先下载到工作目录再识别。',
    importChooseFile: '选择本地文件…',
    importFileChosen: '已选择：{v}',
    importFileUploadHint: '上传文件按自动识别打开；需要指定 member/object 时请改用本地路径导入。',
    importFormat: '格式',
    importFormatsLoading: '正在获取格式列表…',
    importFormatsFallback: '无法获取格式列表，已使用内置列表',
    importRetryFormats: '重试格式列表',
    importAdvanced: '高级选项',
    importTranspose: '矩阵方向',
    importTransposeAuto: '自动判断',
    importTransposeNo: '行=细胞（不转置）',
    importTransposeYes: '行=基因（转置）',
    importObject: 'RData 对象名',
    importObjectPlaceholder: '例如 seu / sce / adata',
    importMember: '压缩包内文件',
    importMemberPlaceholder: '例如 counts/matrix.mtx',
    importSubmit: '开始导入',
    importBusy: '正在导入… {v}',
    importInspect: '预检来源',
    importInspectBusy: '正在预检… {v}',
    importInspectFileHint: '上传的文件无法预检，请改用本地路径或 URL。',
    importInspectResult: '预检结果',
    importInspectFormat: '识别格式',
    importInspectSize: '大小',
    importInspectDims: '维度',
    importNeedSource: '请填写本地路径、URL 或选择文件',
    importCandidates: '检测到多个候选，请选择一个继续',
    importCandidateUploadHint: '上传文件无法指定 member/object：请改用本地路径导入并填写高级项。',
    importSelectCandidate: '选择',
    importClose: '关闭',
    // v3：高级分析
    tabAnalysis: '高级分析',
    analysisTitle: '高级分析',
    analysisSubtitle: '基本 / 组间 / 可视化 / 通路 / 基因集 / 轨迹 / 异质性 / 通讯 / 调控 / 虚拟敲除',
    analysisReload: '刷新清单',
    analysisLoading: '正在获取分析清单…',
    analysisEmpty: '还没有分析清单（Host/引擎未就绪）',
    analysisSelectHint: '从左侧选择一个分析开始',
    analysisNativeShort: '可运行',
    analysisExternalShort: '外部',
    analysisUnavailable: '暂不可用',
    analysisNeeds: '前置条件',
    analysisNeedsHint: '前置提示：{v}',
    analysisRun: '运行分析',
    analysisRunning: '运行中…',
    analysisParams: '参数',
    analysisNoParams: '该分析没有参数',
    analysisResult: '结果',
    analysisMethod: '方法',
    analysisTables: '表格',
    analysisTableRows: '共 {v} 行',
    analysisTableShown: '已显示前 {v} 行',
    analysisFigures: '图',
    analysisFigureNoPreview: '无法预览（文件路径如下）',
    analysisWarnings: '提示',
    analysisNewColumns: '新增 obs 列：{v}（嵌入图着色已可选）',
    analysisExternalReason: '不可用原因',
    analysisExternalAlternative: '可行替代',
    analysisEmptyResult: '没有返回可展示的内容',
    analysisMetrics: '指标',
    analysisGroupOther: '其他',
    analysisNoSelection: '未选择',
    groupBasic: '基本分析',
    groupBetween: '组间对比',
    groupViz: '数据可视化',
    groupPathway: '通路分析',
    groupGeneset: '基因集合分析',
    groupTrajectory: '发育轨迹分析',
    groupHeterogeneity: '肿瘤异质性分析',
    groupCommunication: '细胞通讯分析',
    groupRegulation: '转录调控分析',
    groupKnockout: '虚拟敲除',
    // v4：样本 & 分组
    importPathLines: '本地路径（每行一个，文件=一个样本）',
    importPathLinesPlaceholder: '/Users/you/Downloads/GSM...csv.gz\n/Users/you/Downloads/GSM...csv.gz',
    importMergeDir: '目录合并（每个文件=一个样本）',
    importMultiHint: '多行路径 / 目录合并会合成一个数据集：基因取交集，obs 自动加 sample 列（文件名去扩展名）。',
    importMergedPrefix: '已合并 {v} 个样本',
    importMergedHint: '去「记录」→「样本 & 分组」给样本改名并设置大组（young/old 等）。',
    importGoSamples: '去设置分组',
    samplesTitle: '样本 & 分组',
    samplesSubtitle: '每个源文件=一个样本；改名/设大组后点「保存分组」，会写回 obs 的 sample/group 列。',
    samplesReload: '刷新样本表',
    samplesColSource: '源文件名',
    samplesColCells: '细胞数',
    samplesColName: '样本名',
    samplesColGroup: '大组',
    samplesSave: '保存分组',
    samplesSaving: '保存中…',
    samplesSaved: '已保存',
    samplesUpdatedCells: '影响 {v} 个细胞',
    samplesNewGroup: '新建大组',
    samplesAddGroupPlaceholder: '新大组名，如 young / old',
    samplesAdd: '添加',
    samplesEmpty: '当前数据集没有样本信息（单文件导入时可能没有）。多文件/目录合并导入后会自动出现。',
    samplesNeedDataset: '先载入数据集',
    samplesDirty: '有未保存的修改',
    samplesAutoGroup: '批量正则自动分组',
    samplesAutoPatternPlaceholder: '正则或关键词，如 (GC|LM)、young|old',
    samplesAutoRun: '预览分组',
    samplesAutoPreview: '分组预览（未保存，确认后写入样本表）',
    samplesAutoEmpty: '正则没有匹配到样本',
    samplesApply: '采用并保存',
    samplesCancel: '取消',
    samplesImportMap: '导入样本→大组 CSV',
    samplesImportMapPath: 'CSV 路径（样本名/大组两列）',
    samplesImportMapRun: '读取 CSV',
    samplesImportMapResult: '已读取对照：匹配 {v}',
    samplesJoin: '导入细胞级 metadata',
    samplesJoinPath: '细胞级 CSV 路径（按细胞名 join）',
    samplesJoinKey: '细胞名列（可选）',
    samplesJoinRun: '导入 metadata',
    samplesJoinResult: '已导入列：{v}',
    samplesJoinRate: '匹配率',
    samplesGroupAdderHint: '也可直接在表格最后一列下拉里选现有大组。',
    samplesStatusPrefix: '当前分组列：',
    samplesOneSample: '当前只有 1 个样本；多文件合并导入才会出现样本表。',
    // v5：数据根目录扫描 → 样本表
    scanStepTitle: '第 1 步 · 数据根目录 → 识别样本',
    scanRoot: '数据根目录',
    scanRootRequired: '请先填写数据根目录（GEO 解压后的文件夹）',
    scanRootPlaceholder: '粘贴 GEO 解压后的数据根目录，如 /Users/you/Downloads/GSE246662_RAW',
    scanBrowse: '浏览选择…',
    scanRecursive: '递归扫描所有子目录',
    scanRun: '扫描识别',
    scanRunning: '扫描中…',
    scanErrors: '扫描提示',
    scanEmpty: '还没有扫描结果：填数据根目录后点「扫描识别」。',
    scanSamplesTitle: '识别到的样本（勾选参与分析的样本，可改名和指定大组）',
    scanColSelect: '勾选',
    scanColName: '样本名称',
    scanColRaw: '原始名称',
    scanColFormat: '格式',
    scanColDetail: '详情',
    scanColGroup: '所属分组',
    scanColGeo: 'GEO 说明',
    scanSelectAll: '全选',
    scanClearAll: '清空',
    scanSelectedCount: '已选 {v}',
    scanNewGroup: '新增大组',
    scanNewGroupPlaceholder: '新分组名（如 PT / NT）',
    scanAddGroup: '添加分组',
    scanGroupColor: '颜色',
    scanBatchGroup: '批量分组',
    scanBatchPlaceholder: '目标分组',
    scanApplyBatch: '应用到勾选样本',
    scanGroupPlaceholder: '请选择分组',
    scanDuplicateNames: '样本名重复：{v}（请改成唯一名称）',
    scanConfirm: '确认样本与分组，开始分析',
    scanConfirmBusy: '正在导入…',
    scanGeoTitle: 'GEO 数据库查询（辅助分组）',
    scanGeoPlaceholder: 'GSE 编号，如 GSE246662',
    scanGeoRun: '查询 GEO',
    scanGeoRunning: '查询中…',
    scanGeoMatched: '已按文件名匹配到 {v} 个样本说明',
    scanGeoUnmatched: '未匹配到扫描样本的 GSM',
    scanGeoEmpty: '没有查询到 GEO 样本说明',
    scanAdvancedToggle: '高级：直接打开文件 / URL（单文件 / 多路径 / 目录合并 / 压缩包 / RData…）',
    scanAdvancedHint: '原有导入方式保留在这里，不影响上面的扫描流程。',
    scanBrowseUnavailable: '系统选择框不可用，请手动粘贴数据根目录。',
    scanGroupsEmpty: '还没有分组：先在上面「新增大组」建 young / old 等分组。',
    // v6：降维聚类控制卡
    dimredTitle: '一键降维聚类',
    dimredSub: 'Seurat 6.0 官方流程 + Harmony；结果会覆盖当前 clusters/cell_type/UMAP',
    dimredEngine: '引擎',
    dimredEngineSeurat: 'Seurat 6.0（推荐）',
    dimredEngineBuiltin: '内置引擎',
    dimredHarmony: '开启 Harmony 整合',
    dimredBatchKey: '批次列（batch_key）',
    dimredBatchAuto: '自动（优先 sample 列）',
    dimredBatchHint: '多样本数据建议用 sample 列做批次',
    dimredResolution: '聚类分辨率',
    dimredHint: '参考 20–30 个 cluster；0.05–1.5',
    dimredNPcs: '主成分数 n_pcs',
    dimredNeighbors: 'UMAP 近邻数',
    dimredTopGenes: '高变基因数',
    dimredRun: '一键降维聚类',
    dimredRunBusy: '预设流程运行中…',
    dimredRecluster: '只重聚类',
    dimredReclusterBusy: '重聚类中…',
    dimredWarnings: '运行提示',
    dimredSeuratHint: 'Seurat 6.0 依赖本机 R 库；不可用时会报错，可切内置引擎。',
    // v6：手动注释工作台
    annotTitle: '手动注释工作台',
    annotSubtitle: '选簇 → 点 marker 基因 → 出图 → 填细胞类型 → 保存注释',
    annotClusters: '簇列表',
    annotMarkerTable: '该簇 markers（点击基因）',
    annotGene: '基因',
    annotFeatureplot: 'FeaturePlot（画到嵌入图）',
    annotDotplot: 'DotPlot',
    annotHeatmap: '热图',
    annotSelectedGene: '选中基因',
    annotNoGeneHint: '先在表里点一个基因',
    annotSave: '保存注释',
    annotSaved: '注释已保存',
    annotMappingTable: '注释表（cluster → cell_type）',
    annotCsv: '粘贴注释表 CSV',
    annotCsvPlaceholder: '每行：cluster,细胞类型（如 0,T cell）',
    annotCsvApply: '应用注释表',
    annotPlotBusy: '出图中…',
    annotPlotNote: '提示：后端 dotplot/heatmap 目前按每簇 top-N marker 出图，自定义 genes 参数待后端支持；FeaturePlot 会按基因给嵌入图着色。',
    annotNeedMarkers: '先运行「标记基因」步骤或点「计算标记基因」。',
    annotNeedCluster: '先完成聚类步骤。',
    // v6：分组维护
    samplesGroupsTitle: '已有分组',
    samplesDeleteGroup: '删除',
    samplesGroupUsed: '该分组还有样本在使用：先把这些样本切到其他组或清空分组，再删。',
    samplesGroupDeleted: '已删除分组',
    samplesClearGroup: '— 清空分组 —',
    samplesRemoveSample: '移除',
    samplesColActions: '操作',
    annotCsvEmpty: '没有解析到注释行（每行：cluster,细胞类型）',
    cancel: '取消',
    annotCsvApplied: '已应用 {v} 条注释',
    samplesRemoveHint: '按样本从数据集中移除细胞（会重置聚类/嵌入结果）。',
    samplesRemoved: '已移除 {v} 个细胞；聚类/嵌入已重置，请重新一键降维聚类。',
    // v6：Agent 输入
    agentAskPlaceholder: '想做什么？如：加 harmony、resolution 调到 0.3、找 marker、组间差异、导出报告',
    agentSend: '发送',
    agentParsed: '解析出 {v} 条可执行建议',
    agentRun: '直接运行',
    agentRunning: '运行中…',
    agentResult: '结果',
    agentFallbackTitle: '没听懂这句',
    agentFallbackHint: '可以把这句话发给左下角 Agent；或改写为「加 harmony」「resolution 调到 0.3」「找 marker」「组间差异」「导出报告」。',
    agentNoIntent: '没有解析到可执行意图',
    agentDone: '完成',
    agentFailed: '执行失败',
  }

  var DICT_EN = {
    nav: 'Single-cell',
    title: 'Single-cell workbench',
    subtitle: 'QC · clustering · markers · embeddings',
    emptyTitle: 'No dataset yet',
    emptyHint: 'Load the built-in demo to try it, or open your own .h5ad by absolute path.',
    loadDemo: 'Load demo data',
    openH5ad: 'Open h5ad',
    pathPlaceholder: '/Users/you/data/sample.h5ad',
    refresh: 'Refresh',
    loading: 'Loading…',
    opening: 'Opening…',
    uploading: 'Uploading…',
    openPath: 'Open',
    uploadFile: 'Upload file…',
    engine: 'Engine',
    engineBuiltin: 'Built-in engine',
    engineScanpy: 'scanpy',
    engineNone: 'Python environment not ready',
    engineNotReady: 'Not ready',
    pipeline: 'Pipeline',
    pipelineHint: 'Click a title to expand parameters; finished steps can be re-run.',
    agentPanel: 'Agent',
    agentHint: 'You can also ask the agent in the conversation to drive this workbench.',
    suggestions: 'Next steps',
    copyInstruction: 'Copy instruction',
    copied: 'Copied',
    activity: 'Agent activity',
    tabEmbedding: 'Embedding',
    tabQc: 'QC',
    tabMarkers: 'Markers',
    tabClusters: 'Clusters',
    tabRecord: 'Record',
    tabLog: 'Log',
    colorBy: 'Color by',
    run: 'Run',
    rerun: 'Re-run',
    running: 'Running…',
    done: 'Done',
    pending: 'Pending',
    needPrereq: 'Needs prerequisite',
    generateQc: 'Generate QC plot',
    generateDotplot: 'Generate dot plot',
    cellTypePlaceholder: 'Cell type',
    noData: 'Load data first',
    stepQc: 'QC metrics',
    stepFilter: 'Filter cells',
    stepNormalize: 'Normalize',
    stepHvg: 'HVG',
    stepPca: 'PCA',
    stepNeighbors: 'Neighbors',
    stepCluster: 'Cluster',
    stepEmbed: 'Embedding',
    stepMarkers: 'Markers',
    stepAnnotate: 'Cell types',
    paramMinGenes: 'Min genes',
    paramMaxGenes: 'Max genes',
    paramMinCounts: 'Min UMIs',
    paramMaxPctMt: 'Max mito (%)',
    paramMinCells: 'Min cells per gene',
    paramTargetSum: 'Target sum',
    paramNTopGenes: 'Top genes',
    paramFlavor: 'HVG flavor',
    paramNPcs: 'PCs',
    paramNNeighbors: 'Neighbors',
    paramMetric: 'Metric',
    paramResolution: 'Resolution',
    paramMethod: 'Method',
    paramLayout: 'Layout',
    paramPerplexity: 'Perplexity',
    paramGroupby: 'Group by',
    paramTopN: 'Top n',
    paramMarkerMethod: 'Test',
    flavorSeurat: 'seurat',
    flavorCellRanger: 'cell_ranger',
    metricCosine: 'cosine',
    metricEuclidean: 'euclidean',
    metricManhattan: 'manhattan',
    methodLouvain: 'louvain',
    methodLeiden: 'leiden',
    layoutTsne: 't-SNE',
    layoutPca: 'PCA',
    layoutUmap: 'UMAP',
    groupbyClusters: 'clusters',
    markerTTest: 't-test',
    markerWilcoxon: 'wilcoxon',
    noEmbeddingTitle: 'No embedding yet',
    noEmbeddingHint: 'Run the embedding step to get an interactive scatter here.',
    runEmbed: 'Run embedding',
    embedNeedPca: 'Requires PCA first',
    embedFailed: 'Failed to load embedding',
    embedEmpty: 'No points to show',
    loadingEmbedding: 'Loading embedding…',
    downsampleNotice: 'Downsampled to {v} points',
    resetView: 'Reset view',
    colorCluster: 'Cluster (clusters)',
    colorCellType: 'Cell type (cell_type)',
    colorNGenes: 'Genes (n_genes)',
    colorPctMt: 'Mito (pct_counts_mt)',
    colorGene: 'Gene expression',
    genePlaceholder: 'Gene name, e.g. CD3D',
    geneEmptyHint: 'Type a gene name and press Enter to color cells by expression.',
    geneApply: 'Color',
    legend: 'Legend',
    pointsUnit: 'cells',
    hoverHint: 'Wheel to zoom · drag to pan · double-click to reset · click legend to highlight',
    qcMetrics: 'Metrics',
    qcEmptyHint: 'Pick metrics, then click "Generate QC plot".',
    qcNeedStep: 'Tip: run the QC step in the left pipeline first.',
    qcFigureAlt: 'QC violin plot',
    plotFailed: 'Plot failed',
    markersNeedCluster: 'Run the cluster step before computing markers.',
    markersEmpty: 'No marker results yet',
    loadMarkers: 'Compute markers',
    reloadMarkers: 'Recompute',
    dotplotAlt: 'Marker dot plot',
    tableCluster: 'Cluster',
    tableGene: 'Gene',
    tableLogfc: 'logFC',
    tablePval: 'adj.P',
    tablePctIn: 'pct_in',
    clusterBarTitle: 'Cells per cluster',
    generateClusterBar: 'Generate cluster bar',
    clusterBarFailed: 'Failed to generate cluster bar',
    annotationTitle: 'Cell type annotation',
    annotationHint: 'Edit names and save; they are written back as cell_type.',
    annotationNeedMapping: 'Fill at least one cell type name in the Clusters tab first',
    annotationCountHint: '{v} cell types filled (edit in Clusters tab)',
    saveAnnotation: 'Save annotation',
    annotationSaved: 'Saved',
    cellsUnit: 'cells',
    genesUnit: 'genes',
    logEmpty: 'No activity yet',
    noSuggestions: 'No suggestions',
    activityEmpty: 'No activity',
    progress: 'Progress',
    statusTitle: 'Status',
    colTime: 'Time',
    colStep: 'Step',
    colParams: 'Params',
    colSummary: 'Summary',
    colMs: 'Duration',
    resultOk: 'ok',
    resultFail: 'failed',
    copyFailed: 'Copy failed (clipboard not allowed)',
    requestFailed: 'Request failed',
    pollFailed: 'Failed to refresh state',
    errorTitle: 'Something went wrong',
    retry: 'Retry',
    crashTitle: 'Workbench display error',
    crashHint: 'One block failed to render. Click retry; data and history are unaffected.',
    recipes: 'Recipes',
    recipePick: 'Pick a recipe',
    runRecipe: 'Run recipe',
    recipeRunning: 'Running…',
    recipeStepOf: 'Step {v}',
    recipeDone: 'Recipe finished {v} steps',
    recipeFailed: 'Recipe failed',
    historyTitle: 'Dataset history',
    historyEmpty: 'No history yet; run a step and the provenance trail appears here.',
    auditTitle: 'Self-check',
    runAudit: 'Run self-check',
    auditRunning: 'Checking…',
    auditEmpty: 'No self-check yet',
    auditScore: 'Score {v}',
    auditOk: 'pass',
    auditWarn: 'warn',
    auditFail: 'fail',
    autoannotateTitle: 'Auto-annotation suggestions',
    runAutoannotate: 'Suggest cell types',
    autoannotateRunning: 'Running…',
    autoannotateEmpty: 'No suggestions yet',
    autoannotatePanel: 'Reference panel',
    autoannotateScore: 'score',
    autoannotateMargin: 'margin',
    autoannotateEvidence: 'Evidence',
    autoannotateAlternative: 'Alternative',
    adopt: 'Adopt',
    adoptAll: 'Adopt all',
    adopted: 'Adopted',
    adoptNeedCluster: 'Run the cluster step before auto-annotation.',
    resolutionTitle: 'Resolution scan',
    resolutionPlaceholder: '0.4,0.8,1.2',
    runScan: 'Scan',
    scanRunning: 'Scanning…',
    scanEmpty: 'No scan results yet',
    scanHint: 'Compare cluster counts and ARI across resolutions on the existing kNN graph; the current clustering is unchanged.',
    scanCurrent: 'current',
    colResolution: 'resolution',
    colClusters: 'clusters',
    colAri: 'ARI',
    exportTitle: 'Export artifacts',
    exportProtocol: 'Protocol',
    exportMethods: 'Methods',
    exportReport: 'Report',
    exportScript: 'Script',
    exportH5ad: 'h5ad',
    exportMarkersCsv: 'Markers CSV',
    exportRunning: 'Exporting…',
    exportDone: 'Exported: ',
    copyPath: 'Copy path',
    recordNeedDataset: 'Load a dataset and the provenance record appears here.',
    h5adNeeded: 'Choose a .h5ad / .h5 file',
    pollErrorPrefix: 'Failed to refresh state',
    importData: 'Import data',
    importTitle: 'Import single-cell data',
    importSubtitle: 'h5ad / 10x / loom / matrix text / GEO series / R (.rds/.RData) / archives / URL',
    importSourcePath: 'Local path',
    importSourceUrl: 'Download URL',
    importSourceFile: 'Local file',
    importPathPlaceholder: '/Users/you/data/sample.h5ad or a directory',
    importPathHint: 'File or directory; directories are scanned for 10x triplets. Relative paths resolve against the Host workdir.',
    importUrlPlaceholder: 'https://ftp.ncbi.nlm.nih.gov/geo/series/.../suppl/xxx.txt.gz',
    importUrlHint: 'http(s)/ftp direct links (GEO supplementary files); downloaded to the workdir first.',
    importChooseFile: 'Choose local file…',
    importFileChosen: 'Selected: {v}',
    importFileUploadHint: 'Uploads are opened with auto-detection; use a local path import if you need member/object.',
    importFormat: 'Format',
    importFormatsLoading: 'Loading formats…',
    importFormatsFallback: 'Could not load formats; using the built-in list',
    importRetryFormats: 'Retry formats',
    importAdvanced: 'Advanced options',
    importTranspose: 'Orientation',
    importTransposeAuto: 'Auto-detect',
    importTransposeNo: 'Rows = cells (no transpose)',
    importTransposeYes: 'Rows = genes (transpose)',
    importObject: 'RData object name',
    importObjectPlaceholder: 'e.g. seu / sce / adata',
    importMember: 'Archive member',
    importMemberPlaceholder: 'e.g. counts/matrix.mtx',
    importSubmit: 'Import',
    importBusy: 'Importing… {v}',
    importInspect: 'Inspect source',
    importInspectBusy: 'Inspecting… {v}',
    importInspectFileHint: 'Uploaded files cannot be inspected; use a local path or URL.',
    importInspectResult: 'Inspect result',
    importInspectFormat: 'Format',
    importInspectSize: 'Size',
    importInspectDims: 'Dimensions',
    importNeedSource: 'Provide a local path, a URL, or choose a file',
    importCandidates: 'Multiple candidates found; pick one to continue',
    importCandidateUploadHint: 'Uploads cannot specify member/object: import via a local path and fill the advanced fields.',
    importSelectCandidate: 'Select',
    importClose: 'Close',
    tabAnalysis: 'Advanced',
    analysisTitle: 'Advanced analyses',
    analysisSubtitle: 'Basic / between-group / visualization / pathway / gene sets / trajectory / heterogeneity / communication / regulation / virtual knockout',
    analysisReload: 'Reload',
    analysisLoading: 'Loading analyses…',
    analysisEmpty: 'No analyses available (Host/engine not ready)',
    analysisSelectHint: 'Pick an analysis on the left',
    analysisNativeShort: 'native',
    analysisExternalShort: 'external',
    analysisUnavailable: 'Unavailable',
    analysisNeeds: 'Prerequisites',
    analysisNeedsHint: 'Missing: {v}',
    analysisRun: 'Run analysis',
    analysisRunning: 'Running…',
    analysisParams: 'Parameters',
    analysisNoParams: 'This analysis has no parameters',
    analysisResult: 'Result',
    analysisMethod: 'Method',
    analysisTables: 'Tables',
    analysisTableRows: '{v} rows',
    analysisTableShown: 'showing first {v}',
    analysisFigures: 'Figures',
    analysisFigureNoPreview: 'Preview unavailable (file path below)',
    analysisWarnings: 'Warnings',
    analysisNewColumns: 'New obs columns: {v} (available in embedding color-by)',
    analysisExternalReason: 'Unavailable because',
    analysisExternalAlternative: 'Alternative',
    analysisEmptyResult: 'Nothing to display',
    analysisMetrics: 'Metrics',
    analysisGroupOther: 'Other',
    analysisNoSelection: 'Not selected',
    groupBasic: 'Basic',
    groupBetween: 'Between groups',
    groupViz: 'Visualization',
    groupPathway: 'Pathway',
    groupGeneset: 'Gene sets',
    groupTrajectory: 'Trajectory',
    groupHeterogeneity: 'Tumor heterogeneity',
    groupCommunication: 'Cell communication',
    groupRegulation: 'Regulation',
    groupKnockout: 'Virtual knockout',
    importPathLines: 'Local paths (one per line, each file = one sample)',
    importPathLinesPlaceholder: '/Users/you/Downloads/GSM...csv.gz\n/Users/you/Downloads/GSM...csv.gz',
    importMergeDir: 'Merge directory (each file = one sample)',
    importMultiHint: 'Multiple paths / directory merge produce one dataset: genes intersected, obs gets a sample column (file stem).',
    importMergedPrefix: 'Merged {v} samples',
    importMergedHint: 'Open Record → Samples & groups to rename samples and assign groups (young/old).',
    importGoSamples: 'Set groups',
    samplesTitle: 'Samples & groups',
    samplesSubtitle: 'Each source file = one sample. Renaming / assigning groups then Save writes sample and group columns back to obs.',
    samplesReload: 'Reload',
    samplesColSource: 'Source file',
    samplesColCells: 'Cells',
    samplesColName: 'Sample name',
    samplesColGroup: 'Group',
    samplesSave: 'Save groups',
    samplesSaving: 'Saving…',
    samplesSaved: 'Saved',
    samplesUpdatedCells: '{v} cells updated',
    samplesNewGroup: 'New group',
    samplesAddGroupPlaceholder: 'New group name, e.g. young / old',
    samplesAdd: 'Add',
    samplesEmpty: 'No sample info in this dataset (single-file imports may not have any). Multi-file or directory merges add it automatically.',
    samplesNeedDataset: 'Load a dataset first',
    samplesDirty: 'Unsaved changes',
    samplesAutoGroup: 'Batch auto-group by regex',
    samplesAutoPatternPlaceholder: 'Regex or keywords, e.g. (GC|LM), young|old',
    samplesAutoRun: 'Preview groups',
    samplesAutoPreview: 'Group preview (not saved; apply writes it into the sample table)',
    samplesAutoEmpty: 'The pattern matched no samples',
    samplesApply: 'Apply & save',
    samplesCancel: 'Cancel',
    samplesImportMap: 'Import sample→group CSV',
    samplesImportMapPath: 'CSV path (sample name / group columns)',
    samplesImportMapRun: 'Read CSV',
    samplesImportMapResult: 'Mapping loaded: {v} matched',
    samplesJoin: 'Import cell-level metadata',
    samplesJoinPath: 'Cell-level CSV path (join by cell name)',
    samplesJoinKey: 'Cell name column (optional)',
    samplesJoinRun: 'Import metadata',
    samplesJoinResult: 'Joined columns: {v}',
    samplesJoinRate: 'Match rate',
    samplesGroupAdderHint: 'You can also pick an existing group from the per-row dropdown.',
    samplesStatusPrefix: 'Group column: ',
    samplesOneSample: 'Only one sample; the sample table appears after a multi-file merge import.',
    scanStepTitle: 'Step 1 · Data root → detected samples',
    scanRoot: 'Data root directory',
    scanRootRequired: 'Enter a data root directory first (the unpacked GEO folder)',
    scanRootPlaceholder: 'Paste the GEO data root, e.g. /Users/you/Downloads/GSE246662_RAW',
    scanBrowse: 'Browse…',
    scanRecursive: 'Scan all subdirectories recursively',
    scanRun: 'Scan',
    scanRunning: 'Scanning…',
    scanErrors: 'Scan notes',
    scanEmpty: 'No scan result yet: enter a data root and click Scan.',
    scanSamplesTitle: 'Detected samples (tick the ones to include, rename and assign groups)',
    scanColSelect: 'Include',
    scanColName: 'Sample name',
    scanColRaw: 'Original name',
    scanColFormat: 'Format',
    scanColDetail: 'Detail',
    scanColGroup: 'Group',
    scanColGeo: 'GEO note',
    scanSelectAll: 'Select all',
    scanClearAll: 'Clear',
    scanSelectedCount: '{v} selected',
    scanNewGroup: 'New group',
    scanNewGroupPlaceholder: 'New group name (e.g. PT / NT)',
    scanAddGroup: 'Add group',
    scanGroupColor: 'Color',
    scanBatchGroup: 'Batch group',
    scanBatchPlaceholder: 'Target group',
    scanApplyBatch: 'Apply to selected',
    scanGroupPlaceholder: 'Choose a group',
    scanDuplicateNames: 'Duplicate sample names: {v} (rename to unique names)',
    scanConfirm: 'Confirm samples & groups, start analysis',
    scanConfirmBusy: 'Importing…',
    scanGeoTitle: 'GEO lookup (helps grouping)',
    scanGeoPlaceholder: 'GSE accession, e.g. GSE246662',
    scanGeoRun: 'Look up GEO',
    scanGeoRunning: 'Looking up…',
    scanGeoMatched: 'Matched GEO notes for {v} samples by file name',
    scanGeoUnmatched: 'GSM entries not matched to scanned samples',
    scanGeoEmpty: 'No GEO sample notes returned',
    scanAdvancedToggle: 'Advanced: open a file / URL directly (single file / multiple paths / directory merge / archives / RData…)',
    scanAdvancedHint: 'The previous import flow is kept here; it does not affect the scan flow above.',
    scanBrowseUnavailable: 'The system picker is unavailable; paste the data root manually.',
    scanGroupsEmpty: 'No groups yet: create young / old with "New group" above.',
    dimredTitle: 'One-click embedding & clustering',
    dimredSub: 'Seurat 6.0 official pipeline + Harmony; overwrites clusters/cell_type/UMAP',
    dimredEngine: 'Engine',
    dimredEngineSeurat: 'Seurat 6.0 (recommended)',
    dimredEngineBuiltin: 'Built-in engine',
    dimredHarmony: 'Harmony integration',
    dimredBatchKey: 'Batch column (batch_key)',
    dimredBatchAuto: 'Auto (prefer sample)',
    dimredBatchHint: 'Multi-sample data: use the sample column as batch',
    dimredResolution: 'Resolution',
    dimredHint: 'Target 20–30 clusters; 0.05–1.5',
    dimredNPcs: 'n_pcs',
    dimredNeighbors: 'UMAP neighbors',
    dimredTopGenes: 'n_top_genes',
    dimredRun: 'Run pipeline',
    dimredRunBusy: 'Running recipe…',
    dimredRecluster: 'Re-cluster only',
    dimredReclusterBusy: 'Re-clustering…',
    dimredWarnings: 'Run notes',
    dimredSeuratHint: 'Seurat 6.0 needs the local R library; switch to built-in if unavailable.',
    annotTitle: 'Manual annotation workbench',
    annotSubtitle: 'Pick a cluster → click marker genes → plot → name the cell type → save',
    annotClusters: 'Clusters',
    annotMarkerTable: 'Markers for this cluster (click a gene)',
    annotGene: 'Gene',
    annotFeatureplot: 'FeaturePlot (on embedding)',
    annotDotplot: 'DotPlot',
    annotHeatmap: 'Heatmap',
    annotSelectedGene: 'Selected gene',
    annotNoGeneHint: 'Click a gene in the table first',
    annotSave: 'Save annotation',
    annotSaved: 'Annotation saved',
    annotMappingTable: 'Annotation table (cluster → cell_type)',
    annotCsv: 'Paste annotation CSV',
    annotCsvPlaceholder: 'One per line: cluster,cell type (e.g. 0,T cell)',
    annotCsvApply: 'Apply annotation table',
    annotPlotBusy: 'Plotting…',
    annotPlotNote: 'Note: dotplot/heatmap currently plot top-N markers per cluster; custom genes need backend support. FeaturePlot colors the embedding by gene.',
    annotNeedMarkers: 'Run the markers step or click "Compute markers" first.',
    annotNeedCluster: 'Run the cluster step first.',
    samplesGroupsTitle: 'Groups',
    samplesDeleteGroup: 'Delete',
    samplesGroupUsed: 'This group still has samples: reassign or clear them before deleting.',
    samplesGroupDeleted: 'Group deleted',
    samplesClearGroup: '— clear group —',
    samplesRemoveSample: 'Remove',
    samplesColActions: 'Action',
    annotCsvEmpty: 'No annotation lines parsed (one per line: cluster,cell type)',
    cancel: 'Cancel',
    annotCsvApplied: 'Applied {v} annotations',
    samplesRemoveHint: 'Remove this sample\u2019s cells from the dataset (resets clustering/embedding).',
    samplesRemoved: 'Removed {v} cells; clustering/embedding were reset, re-run dimensionality reduction.',
    agentAskPlaceholder: 'What do you want? e.g. add harmony, set resolution 0.3, find markers, between-group DE, export report',
    agentSend: 'Send',
    agentParsed: '{v} actionable suggestions',
    agentRun: 'Run',
    agentRunning: 'Running…',
    agentResult: 'Result',
    agentFallbackTitle: 'Could not parse this',
    agentFallbackHint: 'Send it to the Agent, or rephrase as "add harmony", "set resolution 0.3", "find markers", "between-group DE", "export report".',
    agentNoIntent: 'No actionable intent found',
    agentDone: 'Done',
    agentFailed: 'Failed',
  }

  // ==========================================================================
  // 基础工具（纯函数）
  // ==========================================================================

  function clamp(value, min, max) {
    var number = Number(value)
    if (!Number.isFinite(number)) return min
    if (number < min) return min
    if (number > max) return max
    return number
  }

  function normalizeDigits(value, fallback) {
    var number = Number(value)
    if (!Number.isFinite(number)) return fallback
    number = Math.floor(number)
    if (number < 0) return 0
    if (number > 6) return 6
    return number
  }

  function trimZeros(text) {
    var value = String(text)
    if (value.indexOf('.') < 0) return value
    value = value.replace(/0+$/, '')
    value = value.replace(/\.$/, '')
    return value === '-0' ? '0' : value
  }

  function groupDigits(text) {
    var raw = String(text)
    var negative = raw.charAt(0) === '-'
    if (negative) raw = raw.slice(1)
    var parts = raw.split('.')
    var integer = parts[0] || '0'
    var grouped = integer.replace(/\B(?=(\d{3})+(?!\d))/g, ',')
    return (negative ? '-' : '') + (parts.length > 1 && parts[1] ? grouped + '.' + parts[1] : grouped)
  }

  /** 数值格式化：≥1 万/1 亿用中文单位，千位分组；非有限值返回「—」。 */
  function formatNumber(value, digits) {
    if (value === null || value === undefined || value === '') return '—'
    var number = Number(value)
    if (!Number.isFinite(number)) return '—'
    var places = normalizeDigits(digits, 1)
    var abs = Math.abs(number)
    if (abs >= 1e8) return trimZeros((number / 1e8).toFixed(places)) + ' 亿'
    if (abs >= 1e4) return trimZeros((number / 1e4).toFixed(places)) + ' 万'
    if (Number.isInteger(number)) return groupDigits(String(number))
    return groupDigits(trimZeros(number.toFixed(places)))
  }

  /** 百分比格式化：入参是 0–1 的小数（pct_in / 占比）。 */
  function formatPct(fraction, digits) {
    if (fraction === null || fraction === undefined || fraction === '') return '—'
    var number = Number(fraction)
    if (!Number.isFinite(number)) return '—'
    return trimZeros((number * 100).toFixed(normalizeDigits(digits, 1))) + '%'
  }

  /** 校正 P 值：极小值用科学计数法，普通值保留 3 位。 */
  function formatPValue(value) {
    var number = Number(value)
    if (!Number.isFinite(number) || number < 0) return '—'
    if (number === 0) return '0'
    if (number < 1e-4) return number.toExponential(1)
    if (number < 0.01) return number.toFixed(4)
    return number.toFixed(3)
  }

  function formatDuration(ms) {
    if (ms === null || ms === undefined || ms === '') return ''
    var number = Number(ms)
    if (!Number.isFinite(number) || number < 0) return ''
    if (number < 1000) return Math.round(number) + ' ms'
    var seconds = number / 1000
    if (seconds < 10) return trimZeros(seconds.toFixed(1)) + ' s'
    if (seconds < 60) return Math.round(seconds) + ' s'
    return trimZeros((seconds / 60).toFixed(1)) + ' min'
  }

  function formatClock(value) {
    if (value === null || value === undefined || value === '') return ''
    var date = value instanceof Date ? value : new Date(String(value))
    if (!date || Number.isNaN(date.getTime())) return ''
    function pad(number) {
      return number < 10 ? '0' + number : String(number)
    }
    return pad(date.getHours()) + ':' + pad(date.getMinutes()) + ':' + pad(date.getSeconds())
  }

  function plainText(value, max) {
    var text = value === null || value === undefined ? '' : String(value)
    text = text.replace(/\s+/g, ' ').trim()
    var limit = Number.isFinite(max) && max > 0 ? Math.floor(max) : 600
    return text.length > limit ? text.slice(0, limit - 1) + '…' : text
  }

  function formatError(error) {
    if (!error) return ''
    if (typeof error === 'string') return plainText(error, 600)
    var message = error.message ? String(error.message) : String(error)
    return plainText(message, 600)
  }

  function isAbortError(error) {
    if (!error) return false
    if (error.name === 'AbortError') return true
    if (error.code === 20) return true
    return /aborted/i.test(String(error.message || ''))
  }

  // ==========================================================================
  // 步骤与参数
  // ==========================================================================

  function stepTitleOf(stepId) {
    return STEP_TITLES_ZH[String(stepId || '')] || String(stepId || '')
  }

  function doneStepSet(doneSteps) {
    var set = {}
    if (!Array.isArray(doneSteps)) return set
    for (var index = 0; index < doneSteps.length; index++) {
      var id = String(doneSteps[index] || '')
      if (STEP_BY_ID[id]) set[id] = true
    }
    return set
  }

  function isStepReady(stepId, doneSteps) {
    var step = STEP_BY_ID[String(stepId || '')]
    if (!step) return false
    var set = doneStepSet(doneSteps)
    for (var index = 0; index < step.requires.length; index++) {
      if (!set[step.requires[index]]) return false
    }
    return true
  }

  /** 当前可运行的步骤（依赖全满足且自己还没完成），按 STEPS 顺序返回。 */
  function nextRunnableSteps(doneSteps) {
    var set = doneStepSet(doneSteps)
    var output = []
    for (var index = 0; index < STEPS.length; index++) {
      var step = STEPS[index]
      if (set[step.id]) continue
      var ready = true
      for (var dep = 0; dep < step.requires.length; dep++) {
        if (!set[step.requires[dep]]) {
          ready = false
          break
        }
      }
      if (ready) output.push(step)
    }
    return output
  }

  function nextStep(doneSteps) {
    var runnable = nextRunnableSteps(doneSteps)
    return runnable.length ? runnable[0] : null
  }

  function stepStatusOf(stepId, doneSteps, runningStepIds) {
    var id = String(stepId || '')
    if (Array.isArray(runningStepIds) && runningStepIds.indexOf(id) >= 0) return 'running'
    if (doneStepSet(doneSteps)[id]) return 'done'
    if (isStepReady(id, doneSteps)) return 'pending'
    return 'blocked'
  }

  function stepStatusKey(status) {
    if (status === 'done') return 'done'
    if (status === 'running') return 'running'
    if (status === 'blocked') return 'needPrereq'
    return 'pending'
  }

  function defaultParams(stepId) {
    var step = STEP_BY_ID[String(stepId || '')]
    var params = {}
    if (!step) return params
    for (var index = 0; index < step.params.length; index++) {
      var param = step.params[index]
      params[param.key] = param.default
    }
    return params
  }

  /**
   * 把 UI 里的原始输入收成 worker 参数。
   * - 数字为空：默认值为 null 的保留 null（如 max_genes），否则省略让引擎用默认值；
   * - 非数字返回可读错误（行内红字用）。
   */
  function collectStepParams(stepId, values) {
    var step = STEP_BY_ID[String(stepId || '')]
    var source = values && typeof values === 'object' ? values : {}
    var params = {}
    if (!step) return { params: params, error: '' }
    for (var index = 0; index < step.params.length; index++) {
      var param = step.params[index]
      var raw = source[param.key]
      var empty = raw === '' || raw === null || raw === undefined
      if (param.kind === 'number') {
        if (empty) {
          if (param.default === null) params[param.key] = null
          continue
        }
        var numeric = Number(raw)
        if (!Number.isFinite(numeric)) return { params: params, error: '参数「' + param.key + '」需要是数字' }
        params[param.key] = numeric
        continue
      }
      if (param.kind === 'select' || param.kind === 'text') {
        if (empty) continue
        params[param.key] = String(raw)
        continue
      }
    }
    return { params: params, error: '' }
  }

  // ==========================================================================
  // 建议 → 复制指令（给 Agent 的中文指令）
  // ==========================================================================

  function formatParamValue(value) {
    if (value === null) return 'null'
    if (value === undefined) return ''
    if (Array.isArray(value)) {
      var parts = []
      for (var index = 0; index < value.length; index++) parts.push(formatParamValue(value[index]))
      return parts.join('/')
    }
    if (typeof value === 'number') return Number.isFinite(value) ? String(value) : 'NaN'
    if (typeof value === 'boolean') return value ? 'true' : 'false'
    if (typeof value === 'object') {
      try {
        return JSON.stringify(value)
      } catch (error) {
        return '[object]'
      }
    }
    return String(value)
  }

  function formatParamList(params) {
    if (!params || typeof params !== 'object' || Array.isArray(params)) return ''
    var keys = Object.keys(params)
    var parts = []
    for (var index = 0; index < keys.length; index++) {
      parts.push(keys[index] + '=' + formatParamValue(params[keys[index]]))
    }
    return parts.join('，')
  }

  /**
   * 一条 suggestion → 可直接发给 Agent 的中文指令。
   * 带 step 的优先用步骤语义（同 CONTRACT §4.2 的例子）；否则退回标题/动作。
   */
  function buildAgentInstruction(suggestion) {
    var item = suggestion && typeof suggestion === 'object' ? suggestion : {}
    var stepId = String(item.step || '')
    if (STEP_BY_ID[stepId]) {
      var paramsText = formatParamList(item.params)
      return '用单细胞工作台把当前数据集跑到' + stepTitleOf(stepId) + (paramsText ? '，' + paramsText : '') + '。'
    }
    var action = String(item.action || item.id || '')
    if (action === 'demo') return '用单细胞工作台载入内置示例数据。'
    if (action === 'open') return '用单细胞工作台打开一个 h5ad 数据集。'
    if (action === 'export') return '用单细胞工作台导出当前数据集。'
    var title = typeof item.title === 'string' && item.title ? item.title : ''
    if (title) return '用单细胞工作台执行「' + title + '」。'
    return '用单细胞工作台继续当前数据集的分析。'
  }

  // ==========================================================================
  // 坐标变换：world(数据坐标) ↔ screen(画布 CSS 像素)
  // screenY = centerY + offsetY - worldY * scale（y 轴向上，符合散点图习惯）
  // ==========================================================================

  function createView() {
    return { scale: 1, offsetX: 0, offsetY: 0 }
  }

  function viewCenter(width, height) {
    var w = Number(width) > 0 ? Number(width) : 1
    var h = Number(height) > 0 ? Number(height) : 1
    return { x: w / 2, y: h / 2 }
  }

  function worldToScreen(point, view, center) {
    var v = view || createView()
    var c = center || { x: 0, y: 0 }
    var p = point || { x: 0, y: 0 }
    var scale = Number(v.scale) || 1
    return {
      x: c.x + (Number(v.offsetX) || 0) + p.x * scale,
      y: c.y + (Number(v.offsetY) || 0) - p.y * scale,
    }
  }

  function screenToWorld(point, view, center) {
    var v = view || createView()
    var c = center || { x: 0, y: 0 }
    var p = point || { x: 0, y: 0 }
    var scale = Number(v.scale) || 1
    return {
      x: (p.x - c.x - (Number(v.offsetX) || 0)) / scale,
      y: (c.y + (Number(v.offsetY) || 0) - p.y) / scale,
    }
  }

  /** 以屏幕上的某个点为锚缩放：该点下的数据坐标保持不动。 */
  function zoomViewAt(view, screenPoint, factor, center, minScale, maxScale) {
    var v = view || createView()
    var anchor = screenPoint || { x: 0, y: 0 }
    var c = center || { x: 0, y: 0 }
    var ratio = Number(factor)
    if (!Number.isFinite(ratio) || ratio <= 0) ratio = 1
    var low = Number.isFinite(minScale) && minScale > 0 ? minScale : 0.01
    var high = Number.isFinite(maxScale) && maxScale > low ? maxScale : 1000
    var current = Number(v.scale) || 1
    var scale = clamp(current * ratio, low, high)
    var world = screenToWorld(anchor, v, c)
    return {
      scale: scale,
      offsetX: anchor.x - c.x - world.x * scale,
      offsetY: anchor.y - c.y + world.y * scale,
    }
  }

  function panView(view, deltaX, deltaY) {
    var v = view || createView()
    return {
      scale: Number(v.scale) || 1,
      offsetX: (Number(v.offsetX) || 0) + (Number(deltaX) || 0),
      offsetY: (Number(v.offsetY) || 0) + (Number(deltaY) || 0),
    }
  }

  function computeBounds(xs, ys) {
    if (!Array.isArray(xs) || !Array.isArray(ys) || !xs.length || !ys.length) return null
    var count = Math.min(xs.length, ys.length)
    var minX = Infinity
    var maxX = -Infinity
    var minY = Infinity
    var maxY = -Infinity
    var found = false
    for (var index = 0; index < count; index++) {
      var x = Number(xs[index])
      var y = Number(ys[index])
      if (!Number.isFinite(x) || !Number.isFinite(y)) continue
      found = true
      if (x < minX) minX = x
      if (x > maxX) maxX = x
      if (y < minY) minY = y
      if (y > maxY) maxY = y
    }
    if (!found) return null
    return { minX: minX, maxX: maxX, minY: minY, maxY: maxY }
  }

  /** 让整幅数据落在 width×height 内（padding 为 CSS 像素），并居中。 */
  function fitView(bounds, width, height, padding) {
    var w = Number(width) > 0 ? Number(width) : 1
    var h = Number(height) > 0 ? Number(height) : 1
    var pad = Number.isFinite(padding) && padding >= 0 ? padding : 24
    if (!bounds) return createView()
    var spanX = bounds.maxX - bounds.minX
    var spanY = bounds.maxY - bounds.minY
    if (!(spanX > 0)) spanX = 1
    if (!(spanY > 0)) spanY = 1
    var availableW = Math.max(1, w - pad * 2)
    var availableH = Math.max(1, h - pad * 2)
    var scale = Math.min(availableW / spanX, availableH / spanY)
    if (!Number.isFinite(scale) || scale <= 0) scale = 1
    var centerX = (bounds.minX + bounds.maxX) / 2
    var centerY = (bounds.minY + bounds.maxY) / 2
    return { scale: scale, offsetX: -centerX * scale, offsetY: centerY * scale }
  }

  /** 命中测试：返回离屏幕点最近的数据点下标（超出 maxDistance 返回 -1）。 */
  function findNearestPoint(xs, ys, view, center, screenPoint, maxDistance) {
    if (!Array.isArray(xs) || !Array.isArray(ys) || !xs.length || !ys.length) return -1
    var count = Math.min(xs.length, ys.length)
    var v = view || createView()
    var c = center || { x: 0, y: 0 }
    var target = screenPoint || { x: 0, y: 0 }
    var limit = Number.isFinite(maxDistance) && maxDistance > 0 ? maxDistance : 12
    var best = -1
    var bestSquared = limit * limit
    for (var index = 0; index < count; index++) {
      var projected = worldToScreen({ x: xs[index], y: ys[index] }, v, c)
      var dx = projected.x - target.x
      var dy = projected.y - target.y
      var squared = dx * dx + dy * dy
      if (squared < bestSquared) {
        bestSquared = squared
        best = index
      }
    }
    return best
  }

  /** 悬停提示文本：分类图给出分组名，连续图给出指标名 + 数值。 */
  function tooltipText(input) {
    var point = input && typeof input === 'object' ? input : {}
    var index = Number(point.index)
    var head = '细胞 #' + (Number.isFinite(index) && index >= 0 ? Math.floor(index) + 1 : '?')
    if (point.continuous) {
      var title = typeof point.title === 'string' && point.title ? point.title : ''
      var valueText = formatNumber(point.value, 2)
      return title ? head + ' · ' + title + ' ' + valueText : head + ' · ' + valueText
    }
    var label = point.label === null || point.label === undefined ? '' : String(point.label)
    return label ? head + ' · ' + label : head
  }

  // ==========================================================================
  // 数据色：连续着色色带插值
  // ==========================================================================

  function hexToRgb(hex) {
    var text = String(hex || '').replace('#', '')
    if (text.length === 3) text = text.charAt(0) + text.charAt(0) + text.charAt(1) + text.charAt(1) + text.charAt(2) + text.charAt(2)
    if (text.length !== 6) return { r: 0, g: 0, b: 0 }
    var value = parseInt(text, 16)
    if (!Number.isFinite(value)) return { r: 0, g: 0, b: 0 }
    return { r: (value >> 16) & 255, g: (value >> 8) & 255, b: value & 255 }
  }

  function lerpColor(from, to, ratio) {
    var a = hexToRgb(from)
    var b = hexToRgb(to)
    var k = clamp(ratio, 0, 1)
    var r = Math.round(a.r + (b.r - a.r) * k)
    var g = Math.round(a.g + (b.g - a.g) * k)
    var blue = Math.round(a.b + (b.b - a.b) * k)
    return 'rgb(' + r + ',' + g + ',' + blue + ')'
  }

  function pickContinuousColor(value, min, max) {
    var number = Number(value)
    if (!Number.isFinite(number)) return CONTINUOUS_STOPS[0]
    var low = Number(min)
    var high = Number(max)
    var t = 0.5
    if (Number.isFinite(low) && Number.isFinite(high) && high > low) t = clamp((number - low) / (high - low), 0, 1)
    var scaled = t * (CONTINUOUS_STOPS.length - 1)
    var index = Math.min(CONTINUOUS_STOPS.length - 2, Math.floor(scaled))
    return lerpColor(CONTINUOUS_STOPS[index], CONTINUOUS_STOPS[index + 1], scaled - index)
  }

  // ==========================================================================
  // 归一化：把 Host 返回的形状收敛成组件能安全消费的形状（缺字段也不崩）
  // ==========================================================================

  function normalizeEngine(engine) {
    var raw = engine && typeof engine === 'object' ? engine : {}
    var mode = raw.mode
    if (mode !== 'scanpy' && mode !== 'builtin' && mode !== 'none' && mode !== 'unknown') {
      mode = raw.ready === true ? 'builtin' : 'unknown'
    }
    return {
      ready: raw.ready === true,
      mode: mode,
      python: typeof raw.python === 'string' ? raw.python : '',
      version: typeof raw.version === 'string' ? raw.version : '',
      packages: raw.packages && typeof raw.packages === 'object' ? raw.packages : {},
      layouts: Array.isArray(raw.layouts) ? raw.layouts.map(String) : [],
      warnings: Array.isArray(raw.warnings) ? raw.warnings.map(String) : [],
      error: typeof raw.error === 'string' ? raw.error : '',
      progress: typeof raw.progress === 'string' ? raw.progress : '',
    }
  }

  function chooseEngine(primary, fallback) {
    var a = primary && typeof primary === 'object' ? primary : null
    var b = fallback && typeof fallback === 'object' ? fallback : null
    function describe(engine) {
      if (!engine) return 0
      if (engine.ready) return 3
      if (engine.mode === 'scanpy' || engine.mode === 'builtin' || engine.mode === 'none') return 2
      return 1
    }
    return describe(b) > describe(a) ? b : (a || b || { ready: false, mode: 'unknown' })
  }

  /** 引擎徽章的展示信息：builtin / scanpy / 未就绪。 */
  function engineBadge(engine) {
    var e = engine && typeof engine === 'object' ? engine : {}
    if (e.mode === 'scanpy') {
      return { key: 'engineScanpy', labelKey: 'engineScanpy', tone: 'ok', detail: e.python ? 'Python ' + e.python : '' }
    }
    if (e.mode === 'builtin') {
      return { key: 'engineBuiltin', labelKey: 'engineBuiltin', tone: 'ok', detail: e.python ? 'Python ' + e.python : '' }
    }
    if (e.mode === 'none' || e.error) {
      return { key: 'engineNone', labelKey: 'engineNotReady', tone: 'bad', detail: e.error || '' }
    }
    return { key: 'engineNone', labelKey: 'engineNotReady', tone: 'idle', detail: e.progress || '' }
  }

  function normalizeStepIds(list) {
    var output = []
    if (!Array.isArray(list)) return output
    for (var index = 0; index < list.length; index++) {
      var id = String(list[index] || '')
      if (STEP_BY_ID[id] && output.indexOf(id) < 0) output.push(id)
    }
    return output
  }

  function normalizeHistory(list) {
    var output = []
    if (!Array.isArray(list)) return output
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object') continue
      output.push({
        at: typeof item.at === 'string' ? item.at : '',
        step: String(item.step || ''),
        params: item.params && typeof item.params === 'object' ? item.params : {},
        ok: item.ok !== false,
        ms: item.ms === null || item.ms === undefined ? null : Number(item.ms),
        summary: typeof item.summary === 'string'
          ? item.summary
          : (item.summary && typeof item.summary.title === 'string' ? item.summary.title : ''),
        nObs: Number(item.n_obs) || 0,
        nVars: Number(item.n_vars) || 0,
        engine: item.engine && typeof item.engine === 'object' ? item.engine : {},
      })
    }
    return output
  }

  function normalizeDataset(dataset) {
    var item = dataset && typeof dataset === 'object' ? dataset : {}
    return {
      id: String(item.id || ''),
      name: typeof item.name === 'string' && item.name ? item.name : String(item.id || ''),
      path: typeof item.path === 'string' ? item.path : '',
      nObs: Number(item.n_obs) || Number(item.nObs) || 0,
      nVars: Number(item.n_vars) || Number(item.nVars) || 0,
      steps: normalizeStepIds(item.steps),
      current: typeof item.current === 'string' ? item.current : null,
      warnings: Array.isArray(item.warnings) ? item.warnings.map(function (warning) { return plainText(warning, 240) }).filter(Boolean) : [],
      obsColumns: normalizeObsColumns(item.obs_columns || item.obsColumns),
      sourceFormat: typeof item.source_format === 'string'
        ? item.source_format
        : (typeof item.sourceFormat === 'string' ? item.sourceFormat : ''),
      notes: Array.isArray(item.notes) ? item.notes.map(String) : [],
      analysesRun: Array.isArray(item.analyses_run)
        ? item.analyses_run.filter(function (entry) { return entry && typeof entry === 'object' })
        : [],
      importState: normalizeImportState(item.import_state || item.importState),
      samples: normalizeSamples(item.samples),
      groups: normalizeStringList(item.groups),
      groupColumn: typeof item.group_column === 'string' && item.group_column ? item.group_column : 'group',
      sampleColumn: typeof item.sample_column === 'string' && item.sample_column ? item.sample_column : 'sample',
      history: normalizeHistory(item.history),
      summaries: item.summaries && typeof item.summaries === 'object' ? item.summaries : {},
    }
  }

  function normalizeSuggestions(list) {
    var output = []
    if (!Array.isArray(list)) return output
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object') continue
      output.push({
        id: String(item.id || item.action || ('sug-' + index)),
        title: typeof item.title === 'string' ? item.title : '',
        detail: typeof item.detail === 'string' ? item.detail : '',
        step: typeof item.step === 'string' ? item.step : '',
        action: typeof item.action === 'string' ? item.action : '',
        params: item.params && typeof item.params === 'object' && !Array.isArray(item.params) ? item.params : {},
      })
    }
    return output
  }

  function normalizeActivity(list) {
    var output = []
    if (!Array.isArray(list)) return output
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object') continue
      output.push({
        at: typeof item.at === 'string' ? item.at : '',
        actor: item.actor === 'agent' ? 'agent' : (item.actor === 'system' ? 'system' : 'user'),
        action: String(item.action || ''),
        dataset: item.dataset === null || item.dataset === undefined ? '' : String(item.dataset),
        detail: typeof item.detail === 'string' ? item.detail : '',
        ok: item.ok !== false,
        ms: item.ms === null || item.ms === undefined ? null : Number(item.ms),
      })
    }
    return output
  }

  function normalizeRecipes(list) {
    var output = []
    if (!Array.isArray(list)) return output
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object' || !item.id) continue
      var steps = []
      if (Array.isArray(item.steps)) {
        for (var stepIndex = 0; stepIndex < item.steps.length; stepIndex++) {
          var step = item.steps[stepIndex]
          if (!step || typeof step !== 'object') continue
          var stepId = String(step.step || step.id || '')
          if (!STEP_BY_ID[stepId]) continue
          steps.push({ step: stepId, params: step.params && typeof step.params === 'object' ? step.params : {} })
        }
      }
      output.push({
        id: String(item.id),
        title: typeof item.title === 'string' ? item.title : String(item.id),
        detail: typeof item.detail === 'string' ? item.detail : '',
        steps: steps,
      })
    }
    return output
  }

  function normalizeRecipeProgress(progress) {
    if (!progress || typeof progress !== 'object') return null
    return {
      recipe: String(progress.recipe || ''),
      title: typeof progress.title === 'string' ? progress.title : '',
      index: Number(progress.index) || 0,
      total: Number(progress.total) || 0,
      step: String(progress.step || ''),
      at: typeof progress.at === 'string' ? progress.at : '',
    }
  }

  function normalizeState(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.state && typeof payload.state === 'object' ? payload.state : payload)
      : {}
    var datasets = Array.isArray(raw.datasets) ? raw.datasets.map(normalizeDataset) : []
    var active = typeof raw.active === 'string' && raw.active
      ? raw.active
      : (datasets.length ? datasets[0].id : null)
    return {
      active: active,
      datasets: datasets,
      engine: normalizeEngine(raw.engine),
      suggestions: normalizeSuggestions(raw.suggestions),
      activity: normalizeActivity(raw.activity),
      recipes: normalizeRecipes(raw.recipes),
      recipeProgress: normalizeRecipeProgress(raw.recipeProgress || raw.recipe_progress),
    }
  }

  function findDataset(state, datasetId) {
    if (!state || !Array.isArray(state.datasets)) return null
    for (var index = 0; index < state.datasets.length; index++) {
      if (state.datasets[index].id === datasetId) return state.datasets[index]
    }
    return null
  }

  function normalizePreview(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.preview && typeof payload.preview === 'object' ? payload.preview : payload)
      : {}
    return {
      genes: Array.isArray(raw.genes) ? raw.genes.map(String) : [],
      obsColumns: normalizeObsColumns(raw.obs_columns),
      nObs: Number(raw.n_obs) || 0,
      nVars: Number(raw.n_vars) || 0,
      steps: normalizeStepIds(raw.steps),
      current: typeof raw.current === 'string' ? raw.current : null,
      notes: Array.isArray(raw.notes) ? raw.notes.map(String) : [],
      importState: normalizeImportState(raw.import_state || raw.importState),
    }
  }

  function normalizeCategory(category, index) {
    var item = category && typeof category === 'object' ? category : {}
    var ordinal = Number(index) || 0
    var id = item.id === null || item.id === undefined ? String(ordinal) : String(item.id)
    return {
      id: id,
      label: item.label === null || item.label === undefined ? id : String(item.label),
      color: typeof item.color === 'string' && item.color ? item.color : PALETTE[ordinal % PALETTE.length],
      n: Number(item.n) || 0,
    }
  }

  function normalizeEmbedding(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.embedding && typeof payload.embedding === 'object' ? payload.embedding : payload)
      : {}
    var xs = Array.isArray(raw.x) ? raw.x : []
    var ys = Array.isArray(raw.y) ? raw.y : []
    var continuous = raw.continuous === true
    var values = Array.isArray(raw.values) ? raw.values : []
    var groupIndex = Array.isArray(raw.group_index) ? raw.group_index : []
    var categories = Array.isArray(raw.categories) ? raw.categories.map(normalizeCategory) : []
    var count = Math.min(xs.length, ys.length)
    if (continuous) count = Math.min(count, values.length)
    else if (groupIndex.length) count = Math.min(count, groupIndex.length)
    var x = []
    var y = []
    var normalizedValues = []
    var groups = []
    for (var index = 0; index < count; index++) {
      var xValue = Number(xs[index])
      var yValue = Number(ys[index])
      x.push(Number.isFinite(xValue) ? xValue : 0)
      y.push(Number.isFinite(yValue) ? yValue : 0)
      if (continuous) {
        var value = Number(values[index])
        normalizedValues.push(Number.isFinite(value) ? value : 0)
      } else {
        var group = Number(groupIndex[index])
        groups.push(Number.isFinite(group) ? group : 0)
      }
    }
    return {
      layout: typeof raw.layout === 'string' ? raw.layout : '',
      nPoints: count,
      downsampled: raw.downsampled === true,
      continuous: continuous,
      legendTitle: typeof raw.legend_title === 'string'
        ? raw.legend_title
        : (typeof raw.legendTitle === 'string' ? raw.legendTitle : ''),
      categories: categories,
      x: x,
      y: y,
      values: normalizedValues,
      groupIndex: groups,
    }
  }

  function normalizeMarkersGroups(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var list = []
    if (Array.isArray(raw.groups)) list = raw.groups
    else if (raw.markers && Array.isArray(raw.markers.groups)) list = raw.markers.groups
    else if (Array.isArray(raw.markers)) list = raw.markers
    var groups = []
    for (var index = 0; index < list.length; index++) {
      var group = list[index] && typeof list[index] === 'object' ? list[index] : {}
      var genes = []
      var rawGenes = Array.isArray(group.genes) ? group.genes : []
      for (var geneIndex = 0; geneIndex < rawGenes.length; geneIndex++) {
        var gene = rawGenes[geneIndex] && typeof rawGenes[geneIndex] === 'object' ? rawGenes[geneIndex] : {}
        genes.push({
          gene: String(gene.gene || gene.name || ''),
          score: Number(gene.score) || 0,
          logfc: Number(gene.logfc === undefined || gene.logfc === null ? gene.logFC : gene.logfc) || 0,
          pvalAdj: Number(gene.pval_adj === undefined || gene.pval_adj === null ? gene.pvalAdj : gene.pval_adj) || 0,
          pctIn: Number(gene.pct_in === undefined || gene.pct_in === null ? gene.pctIn : gene.pct_in) || 0,
          pctOut: Number(gene.pct_out === undefined || gene.pct_out === null ? gene.pctOut : gene.pct_out) || 0,
        })
      }
      groups.push({
        id: group.id === null || group.id === undefined ? String(index) : String(group.id),
        label: group.label === null || group.label === undefined
          ? (group.id === null || group.id === undefined ? String(index) : String(group.id))
          : String(group.label),
        nCells: Number(group.n_cells) || Number(group.nCells) || 0,
        genes: genes,
      })
    }
    return groups
  }

  function normalizeFigure(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.figure && typeof payload.figure === 'object' ? payload.figure : payload)
      : {}
    var url = typeof raw.url === 'string' ? raw.url : ''
    if (!url) return null
    return {
      id: raw.id === null || raw.id === undefined ? '' : String(raw.id),
      url: url,
      width: Number(raw.width) || 0,
      height: Number(raw.height) || 0,
      kind: typeof raw.kind === 'string' ? raw.kind : '',
      file: typeof raw.file === 'string' ? raw.file : '',
    }
  }

  function normalizeAudit(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.audit && typeof payload.audit === 'object' ? payload.audit : payload)
      : {}
    var checks = []
    var list = Array.isArray(raw.checks) ? raw.checks : []
    for (var index = 0; index < list.length; index++) {
      var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
      checks.push({
        id: String(item.id || ('check-' + index)),
        level: item.level === 'fail' ? 'fail' : (item.level === 'warn' ? 'warn' : 'ok'),
        title: typeof item.title === 'string' ? item.title : String(item.id || ''),
        detail: typeof item.detail === 'string' ? item.detail : '',
        hint: typeof item.hint === 'string' ? item.hint : '',
      })
    }
    return {
      score: Number.isFinite(Number(raw.score)) ? Math.round(Number(raw.score)) : 0,
      level: raw.level === 'fail' ? 'fail' : (raw.level === 'warn' ? 'warn' : 'ok'),
      checks: checks,
    }
  }

  function normalizeAutoannotate(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.autoannotate && typeof payload.autoannotate === 'object' ? payload.autoannotate : payload)
      : {}
    var suggestions = []
    var list = Array.isArray(raw.suggestions) ? raw.suggestions : []
    for (var index = 0; index < list.length; index++) {
      var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
      suggestions.push({
        cluster: item.cluster === null || item.cluster === undefined ? String(index) : String(item.cluster),
        label: typeof item.label === 'string' ? item.label : '',
        score: Number(item.score) || 0,
        margin: Number(item.margin) || 0,
        evidence: Array.isArray(item.evidence) ? item.evidence.map(String) : [],
        alternative: typeof item.alternative === 'string' ? item.alternative : '',
      })
    }
    return { panel: typeof raw.panel === 'string' ? raw.panel : '', suggestions: suggestions }
  }

  function normalizeResolutionScan(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var current = null
    if (raw.current && typeof raw.current === 'object') {
      var currentResolution = Number(raw.current.resolution)
      current = {
        resolution: Number.isFinite(currentResolution) ? currentResolution : null,
        nClusters: Number(raw.current.n_clusters) || 0,
      }
    }
    var results = []
    var list = Array.isArray(raw.results) ? raw.results : []
    for (var index = 0; index < list.length; index++) {
      var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
      var resolution = Number(item.resolution)
      var ari = Number(item.ari_vs_current)
      results.push({
        resolution: Number.isFinite(resolution) ? resolution : null,
        nClusters: Number(item.n_clusters) || 0,
        ari: Number.isFinite(ari) ? ari : null,
      })
    }
    return { current: current, results: results }
  }

  function normalizeExportResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    return {
      what: typeof raw.what === 'string' ? raw.what : '',
      file: typeof raw.file === 'string' ? raw.file : '',
      bytes: Number(raw.bytes) || 0,
    }
  }

  // ==========================================================================
  // 记录 tab 用到的纯展示函数
  // ==========================================================================

  /** 一条 history → 展示用字段（时间/步骤/参数摘要/结果摘要/耗时/引擎）。 */
  function formatHistoryEntry(entry) {
    var item = entry && typeof entry === 'object' ? entry : {}
    var stepId = String(item.step || '')
    var summary = ''
    if (typeof item.summary === 'string') summary = item.summary
    else if (item.summary && typeof item.summary === 'object' && typeof item.summary.title === 'string') summary = item.summary.title
    return {
      at: typeof item.at === 'string' ? item.at : '',
      time: formatClock(item.at),
      step: stepId,
      title: stepTitleOf(stepId),
      paramsText: formatParamList(item.params),
      summary: plainText(summary, 200),
      ok: item.ok !== false,
      duration: formatDuration(item.ms),
      engineMode: item.engine && typeof item.engine === 'object' && typeof item.engine.mode === 'string' ? item.engine.mode : '',
      nObs: Number(item.n_obs) || 0,
      nVars: Number(item.n_vars) || 0,
    }
  }

  /** 自动注释建议 → annotate 的 mapping（cluster -> cell type）。 */
  function buildAnnotationMapping(suggestions) {
    var mapping = {}
    if (!Array.isArray(suggestions)) return mapping
    for (var index = 0; index < suggestions.length; index++) {
      var item = suggestions[index] && typeof suggestions[index] === 'object' ? suggestions[index] : {}
      var cluster = item.cluster === null || item.cluster === undefined ? '' : String(item.cluster).trim()
      var label = item.label === null || item.label === undefined ? '' : String(item.label).trim()
      if (!cluster || !label) continue
      if (!Object.prototype.hasOwnProperty.call(mapping, cluster)) mapping[cluster] = label
    }
    return mapping
  }

  /** 「0.4,0.8,1.2」→ 分辨率数组（去重、只留 0.1–5，最多 8 个）。 */
  function parseResolutions(input) {
    var raw = Array.isArray(input) ? input : String(input === null || input === undefined ? '' : input).split(/[,，\s]+/)
    var list = []
    var seen = {}
    for (var index = 0; index < raw.length; index++) {
      var value = Number(raw[index])
      if (!Number.isFinite(value) || value <= 0 || value > 5) continue
      var key = String(value)
      if (seen[key]) continue
      seen[key] = true
      list.push(value)
      if (list.length >= 8) break
    }
    return {
      list: list,
      error: list.length ? '' : '请输入 0.1–5 之间的分辨率，例如 0.4,0.8,1.2',
    }
  }

  function recipeProgressText(progress) {
    if (!progress || typeof progress !== 'object') return ''
    var index = Number(progress.index) || 0
    var total = Number(progress.total) || 0
    var stepId = String(progress.step || '')
    var title = stepTitleOf(stepId)
    var position = total > 0 ? index + '/' + total : String(index)
    return '第 ' + position + ' 步' + (title ? '：' + title : '')
  }

  function activityIcon(actor) {
    if (actor === 'agent') return '🤖'
    if (actor === 'system') return '⚙️'
    return '👤'
  }

  // ==========================================================================
  // v3：全格式导入 + 高级分析层的纯函数
  // ==========================================================================

  /** 高级分析分组（顺序即 UI 顺序；与 Python registry 的 group 字段逐字一致）。 */
  var ANALYSIS_GROUPS = [
    '基本分析', '组间对比', '数据可视化', '通路分析', '基因集合分析',
    '发育轨迹分析', '肿瘤异质性分析', '细胞通讯分析', '转录调控分析', '虚拟敲除',
  ]

  /** 分组名 → 字典 key（分组名本身来自服务端 registry，属于数据；UI 优先走翻译）。 */
  var ANALYSIS_GROUP_LABEL_KEYS = {
    '基本分析': 'groupBasic',
    '组间对比': 'groupBetween',
    '数据可视化': 'groupViz',
    '通路分析': 'groupPathway',
    '基因集合分析': 'groupGeneset',
    '发育轨迹分析': 'groupTrajectory',
    '肿瘤异质性分析': 'groupHeterogeneity',
    '细胞通讯分析': 'groupCommunication',
    '转录调控分析': 'groupRegulation',
    '虚拟敲除': 'groupKnockout',
  }

  var ANALYSIS_GROUP_OTHER = '其他'

  /** 导入文件选择器 accept（覆盖 CONTRACT §8.2 全部格式）。 */
  var IMPORT_ACCEPT = '.h5ad,.h5ad.gz,.h5,.loom,.loom.gz,.mtx,.mtx.gz,.csv,.csv.gz,.tsv,.tsv.gz,.txt,.txt.gz,.gz,.zip,.tar,.tar.gz,.tgz,.tar.bz2,.rds,.rdata,.rda,.RDS,.RData'

  /** GET /formats 不可用时的兜底清单：保证导入面板始终能用（服务端优先）。 */
  var FALLBACK_FORMATS = [
    { id: 'auto', title: '自动识别', extensions: [], available: true, note: '按扩展名与内容自动判断格式' },
    { id: 'h5ad', title: 'h5ad (AnnData)', extensions: ['.h5ad', '.h5ad.gz'], available: true, note: '' },
    { id: '10x_h5', title: '10x CellRanger .h5', extensions: ['.h5'], available: true, note: '' },
    { id: '10x_mtx', title: '10x MTX 三件套', extensions: ['.mtx', '.mtx.gz'], available: true, note: '目录或单文件（自动找 sidecar）' },
    { id: 'loom', title: 'loom', extensions: ['.loom'], available: true, note: '' },
    { id: 'matrix_text', title: '矩阵文本 (csv/tsv/txt)', extensions: ['.csv', '.tsv', '.txt', '.gz'], available: true, note: '' },
    { id: 'geo_series', title: 'GEO series matrix', extensions: ['.txt.gz'], available: true, note: '' },
    { id: 'rds', title: 'R .rds (Seurat/SCE)', extensions: ['.rds', '.RDS'], available: true, note: '需要本机 Rscript' },
    { id: 'rdata', title: 'R .RData/.rda', extensions: ['.RData', '.rdata', '.rda'], available: true, note: '需要本机 Rscript；多对象时选 object' },
    { id: 'archive', title: '压缩包 (zip/tar.gz)', extensions: ['.zip', '.tar', '.tar.gz', '.tgz', '.tar.bz2'], available: true, note: '解压后自动识别，可用 member 指定' },
  ]

  function formatBytes(bytes) {
    if (bytes === null || bytes === undefined || bytes === '') return '—'
    var number = Number(bytes)
    if (!Number.isFinite(number) || number < 0) return '—'
    if (number < 1024) return Math.round(number) + ' B'
    var units = ['KB', 'MB', 'GB', 'TB']
    var value = number / 1024
    var unitIndex = 0
    while (value >= 1024 && unitIndex < units.length - 1) {
      value /= 1024
      unitIndex++
    }
    return trimZeros(value.toFixed(1)) + ' ' + units[unitIndex]
  }

  function normalizeFormats(payload) {
    var raw = payload && typeof payload === 'object'
      ? (Array.isArray(payload.formats) ? payload.formats : (Array.isArray(payload) ? payload : []))
      : []
    var output = []
    for (var index = 0; index < raw.length; index++) {
      var item = raw[index]
      if (!item || typeof item !== 'object' || !item.id) continue
      output.push({
        id: String(item.id),
        title: typeof item.title === 'string' && item.title ? item.title : String(item.id),
        extensions: Array.isArray(item.extensions) ? item.extensions.map(String) : [],
        available: item.available !== false,
        note: typeof item.note === 'string' ? item.note : '',
      })
    }
    if (!output.length) {
      for (var fallbackIndex = 0; fallbackIndex < FALLBACK_FORMATS.length; fallbackIndex++) {
        var fallback = FALLBACK_FORMATS[fallbackIndex]
        output.push({
          id: fallback.id,
          title: fallback.title,
          extensions: fallback.extensions.slice(),
          available: fallback.available !== false,
          note: fallback.note || '',
        })
      }
    }
    return output
  }

  function normalizeObsColumns(columns) {
    var output = []
    if (!Array.isArray(columns)) return output
    for (var index = 0; index < columns.length; index++) {
      var column = columns[index]
      if (!column || typeof column !== 'object') continue
      var key = String(column.key || '')
      if (!key) continue
      output.push({ key: key, kind: column.kind === 'numeric' ? 'numeric' : 'categorical' })
    }
    return output
  }

  function normalizeImportInspect(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.inspect && typeof payload.inspect === 'object' ? payload.inspect : payload)
      : {}
    return {
      format: typeof raw.format === 'string' ? raw.format : '',
      size: Number(raw.size) || 0,
      nObs: Number(raw.n_obs) || 0,
      nVars: Number(raw.n_vars) || 0,
      notes: Array.isArray(raw.notes) ? raw.notes.map(String) : [],
      candidates: formatCandidates(raw).items,
    }
  }

  function candidateString(item, keys) {
    for (var index = 0; index < keys.length; index++) {
      var value = item[keys[index]]
      if (typeof value === 'string' && value.trim()) return value.trim()
    }
    return ''
  }

  function normalizeCandidate(candidate, index) {
    var ordinal = Number(index) || 0
    var item = candidate && typeof candidate === 'object' ? candidate : {}
    if (typeof candidate === 'string') item = { path: candidate }
    else if (typeof candidate === 'number') item = { object: String(candidate) }
    var member = candidateString(item, ['member', 'relative_path', 'relativePath', 'archive_path', 'archivePath', 'file'])
    var objectName = candidateString(item, ['object', 'object_name', 'objectName', 'variable'])
    if (!objectName) {
      var weakName = candidateString(item, ['name'])
      if (weakName) {
        var looksLikeFile = /\.[A-Za-z0-9]{1,10}$/.test(weakName)
        var rClass = typeof item.class === 'string' ? item.class : ''
        var formatText = typeof item.format === 'string' ? item.format : ''
        if (rClass || /^(rds|rdata)$/i.test(formatText) || !looksLikeFile) objectName = weakName
      }
    }
    if (!member && typeof item.path === 'string' && item.path.trim()) member = item.path.trim()
    var dimensionText = ''
    if (Array.isArray(item.dim) && item.dim.length >= 2) {
      dimensionText = formatNumber(item.dim[0]) + '×' + formatNumber(item.dim[1])
    } else if (Number(item.n_obs) > 0 && Number(item.n_vars) > 0) {
      dimensionText = formatNumber(item.n_obs) + '×' + formatNumber(item.n_vars)
    }
    var detailParts = []
    var classText = candidateString(item, ['class', 'kind', 'type', 'format'])
    if (classText) detailParts.push(classText)
    if (dimensionText) detailParts.push(dimensionText)
    var size = Number(item.size) || Number(item.bytes) || 0
    if (size > 0) detailParts.push(formatBytes(size))
    var label = candidateString(item, ['label', 'title'])
    if (!label) label = member || objectName || ('候选 ' + (ordinal + 1))
    var patch = {}
    if (member) patch.member = member
    if (objectName) patch.object = objectName
    return {
      key: 'cand-' + ordinal + '-' + label,
      label: label,
      detail: detailParts.join(' · '),
      patch: patch,
      raw: item,
    }
  }

  /** ambiguous 响应（或 candidates 数组）→ 可直接渲染的候选列表。 */
  function formatCandidates(input) {
    var raw = input
    if (raw && typeof raw === 'object' && !Array.isArray(raw)) {
      if (Array.isArray(raw.candidates)) raw = raw.candidates
      else if (raw.data && typeof raw.data === 'object' && Array.isArray(raw.data.candidates)) raw = raw.data.candidates
      else if (raw.error && typeof raw.error === 'object' && Array.isArray(raw.error.candidates)) raw = raw.error.candidates
      else raw = []
    }
    if (!Array.isArray(raw)) raw = []
    var items = []
    for (var index = 0; index < raw.length; index++) items.push(normalizeCandidate(raw[index], index))
    return { items: items, count: items.length }
  }

  /** 导入表单 → POST /dataset/open | /import/inspect 的 body（路径/URL 二选一）。 */
  function buildImportBody(options) {
    var source = options && typeof options === 'object' ? options : {}
    var body = { actor: 'user' }
    var filePath = typeof source.path === 'string' ? source.path.trim() : ''
    var url = typeof source.url === 'string' ? source.url.trim() : ''
    var paths = Array.isArray(source.paths)
      ? source.paths.map(function (item) { return String(item || '').trim() }).filter(Boolean)
      : []
    if (paths.length > 1) body.paths = paths
    else if (filePath) body.path = filePath
    else if (paths.length === 1) body.path = paths[0]
    if (url) body.url = url
    if (source.mergeDir === true || source.merge_dir === true) body.merge_dir = true
    if (Array.isArray(source.samples) && source.samples.length) {
      body.samples = source.samples.map(function (item) {
        var entry = {
          path: String(item && item.path ? item.path : ''),
          name: item && item.name ? String(item.name) : undefined,
          group: item && item.group ? String(item.group) : undefined,
          trio_prefix: item && item.trio_prefix !== undefined && item.trio_prefix !== null ? String(item.trio_prefix) : undefined,
          format: item && item.format ? String(item.format) : undefined,
        }
        return entry
      }).filter(function (item) { return item.path })
    }
    if (typeof source.name === 'string' && source.name.trim()) body.name = source.name.trim()
    var format = typeof source.format === 'string' ? source.format.trim() : ''
    if (format && format !== 'auto') {
      // Host 的 HTTP/工具两条路径分别读 format 与 data_format；两个都发最稳。
      body.format = format
      body.data_format = format
    }
    var member = typeof source.member === 'string' ? source.member.trim() : ''
    if (member) body.member = member
    var objectName = typeof source.object === 'string' ? source.object.trim() : ''
    if (objectName) body.object = objectName
    if (source.transpose === true || source.transpose === 'true') body.transpose = true
    else if (source.transpose === false || source.transpose === 'false') body.transpose = false
    return body
  }

  function normalizeAnalysisParam(param) {
    var item = param && typeof param === 'object' ? param : {}
    var type = item.type
    if (type !== 'number' && type !== 'string' && type !== 'text' && type !== 'select' &&
        type !== 'boolean' && type !== 'column' && type !== 'file') {
      type = 'string'
    }
    var options = []
    if (Array.isArray(item.options)) {
      for (var index = 0; index < item.options.length; index++) {
        var option = item.options[index]
        if (option && typeof option === 'object') {
          var optionValue = option.value === undefined || option.value === null
            ? (option.label === undefined ? '' : String(option.label))
            : String(option.value)
          options.push({
            value: optionValue,
            label: option.label === undefined || option.label === null ? optionValue : String(option.label),
          })
        } else if (option !== null && option !== undefined) {
          options.push({ value: String(option), label: String(option) })
        }
      }
    }
    return {
      key: String(item.key || ''),
      label: typeof item.label === 'string' && item.label ? item.label : String(item.key || ''),
      type: type,
      default: item.default,
      options: options,
      hint: typeof item.hint === 'string' ? item.hint : '',
      min: Number.isFinite(Number(item.min)) ? Number(item.min) : null,
      max: Number.isFinite(Number(item.max)) ? Number(item.max) : null,
      step: Number.isFinite(Number(item.step)) ? Number(item.step) : null,
    }
  }

  function normalizeAnalysisParams(params) {
    var output = []
    if (!Array.isArray(params)) return output
    for (var index = 0; index < params.length; index++) {
      var descriptor = normalizeAnalysisParam(params[index])
      if (!descriptor.key) continue
      output.push(descriptor)
    }
    return output
  }

  /** 合并多个 obs 列清单（preview 带 kind 优先；analyses/serialized 里只有 key）。 */
  function mergeObsColumns() {
    var output = []
    var seen = {}
    for (var argIndex = 0; argIndex < arguments.length; argIndex++) {
      var list = normalizeObsColumns(arguments[argIndex])
      for (var index = 0; index < list.length; index++) {
        var column = list[index]
        if (seen[column.key]) continue
        seen[column.key] = true
        output.push(column)
      }
    }
    return output
  }

  /** GET /analyses 的完整载荷：目录项 + 分组顺序 + obs 列清单 + 诊断信息。 */
  function normalizeAnalysesPayload(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    return {
      dataset: raw.dataset === null || raw.dataset === undefined ? null : String(raw.dataset),
      groups: Array.isArray(raw.groups) ? raw.groups.map(String) : [],
      columns: normalizeObsColumns(raw.columns),
      groupKeys: Array.isArray(raw.group_keys) ? raw.group_keys.map(String) : [],
      diagnostics: raw.diagnostics && typeof raw.diagnostics === 'object' && !Array.isArray(raw.diagnostics) ? raw.diagnostics : {},
      analyses: normalizeAnalyses(payload),
    }
  }

  function normalizeAnalyses(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var list = Array.isArray(raw.analyses)
      ? raw.analyses
      : (Array.isArray(raw.items) ? raw.items : (Array.isArray(payload) ? payload : []))
    var output = []
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object' || !item.id) continue
      var native = item.native !== false
      output.push({
        id: String(item.id),
        title: typeof item.title === 'string' && item.title ? item.title : String(item.id),
        group: typeof item.group === 'string' && item.group ? item.group : ANALYSIS_GROUP_OTHER,
        method: typeof item.method === 'string' ? item.method : '',
        native: native,
        external: !native,
        available: native && item.available !== false,
        reason: typeof item.reason === 'string' ? item.reason : '',
        alternative: typeof item.alternative === 'string' ? item.alternative : '',
        needs: Array.isArray(item.needs) ? item.needs.map(String) : [],
        params: normalizeAnalysisParams(item.params),
      })
    }
    return output
  }

  /** 按 ANALYSIS_GROUPS 顺序分组；未登记的分组按出现顺序排在最后。 */
  function groupAnalyses(list) {
    var items = Array.isArray(list) ? list : []
    var buckets = {}
    var order = []
    for (var index = 0; index < items.length; index++) {
      var item = items[index]
      if (!item || typeof item !== 'object') continue
      var group = typeof item.group === 'string' && item.group ? item.group : ANALYSIS_GROUP_OTHER
      if (!buckets[group]) {
        buckets[group] = []
        order.push(group)
      }
      buckets[group].push(item)
    }
    var groups = []
    var used = {}
    for (var knownIndex = 0; knownIndex < ANALYSIS_GROUPS.length; knownIndex++) {
      var known = ANALYSIS_GROUPS[knownIndex]
      if (buckets[known] && buckets[known].length) {
        groups.push({ group: known, labelKey: ANALYSIS_GROUP_LABEL_KEYS[known] || '', items: buckets[known] })
        used[known] = true
      }
    }
    for (var orderIndex = 0; orderIndex < order.length; orderIndex++) {
      var groupName = order[orderIndex]
      if (used[groupName]) continue
      groups.push({ group: groupName, labelKey: ANALYSIS_GROUP_LABEL_KEYS[groupName] || '', items: buckets[groupName] })
    }
    return groups
  }

  function paramOptionValues(param) {
    return normalizeAnalysisParam(param).options
  }

  function paramDefaultValue(param) {
    var descriptor = normalizeAnalysisParam(param)
    if (descriptor.type === 'number') {
      return Number.isFinite(Number(descriptor.default)) && descriptor.default !== null && descriptor.default !== ''
        ? Number(descriptor.default)
        : ''
    }
    if (descriptor.type === 'boolean') return descriptor.default === true
    if (descriptor.type === 'select') {
      var wanted = descriptor.default === undefined || descriptor.default === null ? '' : String(descriptor.default)
      if (wanted) {
        for (var index = 0; index < descriptor.options.length; index++) {
          if (descriptor.options[index].value === wanted) return wanted
        }
      }
      return descriptor.options.length ? descriptor.options[0].value : ''
    }
    if (descriptor.default === undefined || descriptor.default === null) return ''
    return String(descriptor.default)
  }

  function paramDefaults(params, obsColumns, preferredColumn) {
    var descriptors = normalizeAnalysisParams(params)
    var output = {}
    var columnKeys = {}
    var columns = normalizeObsColumns(obsColumns)
    for (var columnIndex = 0; columnIndex < columns.length; columnIndex++) columnKeys[columns[columnIndex].key] = true
    var preferred = typeof preferredColumn === 'string' && preferredColumn ? preferredColumn : ''
    for (var index = 0; index < descriptors.length; index++) {
      var descriptor = descriptors[index]
      var value = paramDefaultValue(descriptor)
      // v4：column 参数默认值为空时，优先落到数据集的分组列（group_column）
      if (descriptor.type === 'column' && preferred && !value && columnKeys[preferred]) {
        value = preferred
      }
      output[descriptor.key] = value
    }
    return output
  }

  /**
   * 分析参数表单 → POST /analysis 的 params。
   * 空值省略让引擎用默认值；数字非法返回可读错误（行内红字）。
   */
  function collectAnalysisParams(params, values) {
    var descriptors = normalizeAnalysisParams(params)
    var source = values && typeof values === 'object' ? values : {}
    var output = {}
    for (var index = 0; index < descriptors.length; index++) {
      var descriptor = descriptors[index]
      var raw = source[descriptor.key]
      if (descriptor.type === 'number') {
        if (raw === '' || raw === null || raw === undefined) continue
        var number = Number(raw)
        if (!Number.isFinite(number)) {
          return { params: output, error: '参数「' + descriptor.label + '」需要是数字' }
        }
        output[descriptor.key] = number
        continue
      }
      if (descriptor.type === 'boolean') {
        output[descriptor.key] = raw === true || raw === 'true' || raw === 1 || raw === '1'
        continue
      }
      if (raw === '' || raw === null || raw === undefined) continue
      output[descriptor.key] = String(raw)
    }
    return { params: output, error: '' }
  }

  /**
   * needs → 给用户的前置提示。
   * dataset 形状：{steps:[...], obsColumns:[{key,kind}]}（兼容 obs_columns）。
   */
  function needsHint(item, dataset) {
    var analysis = item && typeof item === 'object' ? item : {}
    var needs = Array.isArray(analysis.needs) ? analysis.needs : []
    if (!needs.length) return ''
    var source = dataset && typeof dataset === 'object' ? dataset : {}
    var steps = Array.isArray(source.steps) ? source.steps : []
    var columns = Array.isArray(source.obsColumns)
      ? source.obsColumns
      : (Array.isArray(source.obs_columns) ? source.obs_columns : [])
    var columnKeys = {}
    for (var columnIndex = 0; columnIndex < columns.length; columnIndex++) {
      var column = columns[columnIndex]
      if (column && column.key) columnKeys[String(column.key)] = true
    }
    var hints = []
    for (var index = 0; index < needs.length; index++) {
      var need = String(needs[index] || '')
      if (!need) continue
      if (need === 'clusters') {
        if (steps.indexOf('cluster') < 0) hints.push('需要先运行聚类步骤')
        continue
      }
      if (need.indexOf('groups:') === 0) {
        var groupColumn = need.slice('groups:'.length)
        if (!columnKeys[groupColumn]) hints.push('需要 obs 列「' + groupColumn + '」（导入时提供或先做注释/分组）')
        continue
      }
      if (need === 'gene_sets') {
        hints.push('需要基因集（默认使用内置基因集，也可选 .gmt 文件）')
        continue
      }
      if (need === 'spliced') {
        hints.push('需要 spliced/unspliced 层（当前数据没有，相关分析只能登记为外部工具）')
        continue
      }
      if (need === 'gene_order') {
        hints.push('需要基因位置文件（可选；没有时用不带染色体排序的代理分数）')
        continue
      }
      hints.push('需要 ' + need)
    }
    return hints.join('；')
  }

  function formatTableCell(value, numeric) {
    if (value === null || value === undefined || value === '') return '—'
    if (typeof value === 'number') {
      if (!Number.isFinite(value)) return '—'
      return numeric === false ? trimZeros(String(value)) : formatNumber(value, 3)
    }
    if (typeof value === 'boolean') return value ? 'true' : 'false'
    if (Array.isArray(value)) return value.join(', ')
    if (typeof value === 'object') {
      try {
        return plainText(JSON.stringify(value), 160)
      } catch (error) {
        return '[object]'
      }
    }
    return String(value)
  }

  /** 把 analysis 返回的 table 变成可渲染的 {columns, rows(cells), total, shown, truncated}。 */
  function formatTableRows(table, maxRows) {
    var raw = table && typeof table === 'object' ? table : {}
    var columns = []
    var rawColumns = Array.isArray(raw.columns) ? raw.columns : []
    for (var columnIndex = 0; columnIndex < rawColumns.length; columnIndex++) {
      var column = rawColumns[columnIndex]
      if (!column || typeof column !== 'object') continue
      var key = String(column.key || '')
      if (!key) continue
      columns.push({
        key: key,
        label: typeof column.label === 'string' && column.label ? column.label : key,
        numeric: column.type === 'number',
      })
    }
    var rows = Array.isArray(raw.rows) ? raw.rows : []
    var limit = Number.isFinite(maxRows) && maxRows > 0 ? Math.floor(maxRows) : 0
    var shownRows = limit > 0 ? rows.slice(0, limit) : rows
    var cells = []
    for (var rowIndex = 0; rowIndex < shownRows.length; rowIndex++) {
      var row = shownRows[rowIndex]
      var record = row && typeof row === 'object' && !Array.isArray(row)
        ? row
        : (Array.isArray(row) ? row : {})
      var rowCells = []
      for (var cellIndex = 0; cellIndex < columns.length; cellIndex++) {
        var descriptor = columns[cellIndex]
        var value = Array.isArray(row) ? row[cellIndex] : record[descriptor.key]
        rowCells.push({
          key: descriptor.key,
          text: formatTableCell(value, descriptor.numeric),
          numeric: descriptor.numeric,
        })
      }
      cells.push({ key: 'row-' + rowIndex, cells: rowCells })
    }
    return {
      id: typeof raw.id === 'string' ? raw.id : '',
      title: typeof raw.title === 'string' ? raw.title : '',
      note: typeof raw.note === 'string' ? raw.note : '',
      columns: columns,
      rows: cells,
      total: rows.length,
      shown: cells.length,
      truncated: limit > 0 && rows.length > limit,
    }
  }

  function normalizeAnalysisTable(table) {
    var raw = table && typeof table === 'object' ? table : {}
    var columns = []
    var rawColumns = Array.isArray(raw.columns) ? raw.columns : []
    for (var index = 0; index < rawColumns.length; index++) {
      var column = rawColumns[index]
      if (!column || typeof column !== 'object') continue
      columns.push({
        key: String(column.key || ''),
        label: typeof column.label === 'string' && column.label ? column.label : String(column.key || ''),
        type: column.type === 'number' ? 'number' : 'string',
      })
    }
    return {
      id: typeof raw.id === 'string' ? raw.id : '',
      title: typeof raw.title === 'string' ? raw.title : '',
      note: typeof raw.note === 'string' ? raw.note : '',
      columns: columns,
      rows: Array.isArray(raw.rows) ? raw.rows : [],
    }
  }

  function normalizeAnalysisFigure(figure) {
    var raw = figure && typeof figure === 'object' ? figure : {}
    var url = typeof raw.url === 'string' ? raw.url : ''
    var file = typeof raw.file === 'string' ? raw.file : ''
    if (!url && !file) return null
    return {
      id: raw.id === null || raw.id === undefined ? '' : String(raw.id),
      title: typeof raw.title === 'string' ? raw.title : '',
      url: url,
      file: file,
      width: Number(raw.width) || 0,
      height: Number(raw.height) || 0,
    }
  }

  function normalizeAnalysisResult(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.analysis && typeof payload.analysis === 'object'
          ? payload.analysis
          : (payload.result && typeof payload.result === 'object' ? payload.result : payload))
      : {}
    var summary = raw.summary && typeof raw.summary === 'object' ? raw.summary : {}
    var tables = []
    if (Array.isArray(raw.tables)) {
      for (var tableIndex = 0; tableIndex < raw.tables.length; tableIndex++) {
        tables.push(normalizeAnalysisTable(raw.tables[tableIndex]))
      }
    }
    var figures = []
    if (Array.isArray(raw.figures)) {
      for (var figureIndex = 0; figureIndex < raw.figures.length; figureIndex++) {
        var figure = normalizeAnalysisFigure(raw.figures[figureIndex])
        if (figure) figures.push(figure)
      }
    }
    return {
      id: typeof raw.id === 'string' ? raw.id : '',
      title: typeof raw.title === 'string' ? raw.title : '',
      method: typeof raw.method === 'string' ? raw.method : '',
      summary: {
        bullets: Array.isArray(summary.bullets) ? summary.bullets.map(String) : [],
        metrics: summary.metrics && typeof summary.metrics === 'object' && !Array.isArray(summary.metrics)
          ? summary.metrics
          : {},
      },
      tables: tables,
      figures: figures,
      obsColumns: normalizeObsColumns(raw.obs_columns),
      obsMerged: Array.isArray(raw.obs_merged) ? raw.obs_merged.map(String) : [],
      warnings: Array.isArray(raw.warnings) ? raw.warnings.map(String) : [],
    }
  }

  /** 嵌入图 color-by 选项：固定 4 项 + 数据集的 obs 列（分析新增的列自动出现）+ 基因。 */
  var FIXED_COLOR_OPTIONS = [
    { value: 'cluster:clusters', labelKey: 'colorCluster' },
    { value: 'obs:cell_type', labelKey: 'colorCellType' },
    { value: 'obs:n_genes', labelKey: 'colorNGenes' },
    { value: 'obs:pct_counts_mt', labelKey: 'colorPctMt' },
  ]

  function buildColorOptions(obsColumns) {
    var options = []
    var seen = {}
    for (var fixedIndex = 0; fixedIndex < FIXED_COLOR_OPTIONS.length; fixedIndex++) {
      var fixed = FIXED_COLOR_OPTIONS[fixedIndex]
      options.push({ value: fixed.value, labelKey: fixed.labelKey })
    }
    seen.clusters = true
    seen.cell_type = true
    seen.n_genes = true
    seen.pct_counts_mt = true
    var columns = normalizeObsColumns(obsColumns)
    var extras = []
    for (var index = 0; index < columns.length; index++) {
      var column = columns[index]
      if (seen[column.key]) continue
      seen[column.key] = true
      extras.push({
        value: 'obs:' + column.key,
        label: column.key + (column.kind === 'numeric' ? ' · numeric' : ' · categorical'),
      })
    }
    extras.sort(function (left, right) { return left.value < right.value ? -1 : (left.value > right.value ? 1 : 0) })
    for (var extraIndex = 0; extraIndex < extras.length; extraIndex++) options.push(extras[extraIndex])
    options.push({ value: 'gene', labelKey: 'colorGene' })
    return options
  }

  // ==========================================================================
  // v4：样本 & 分组（GEO 每个 GSM 一个文件 → 合并成多样本数据集）
  // ==========================================================================

  /** 本地路径多行输入 → 路径列表（去空行、去首尾空白、支持逗号/中文逗号分隔）。 */
  function parsePathLines(text) {
    var raw = text === null || text === undefined ? '' : String(text)
    var lines = raw.split(/\r?\n|,|，/)
    var paths = []
    var seen = {}
    for (var index = 0; index < lines.length; index++) {
      var value = String(lines[index] || '').trim()
      if (!value || seen[value]) continue
      seen[value] = true
      paths.push(value)
    }
    return paths
  }

  function normalizeSamples(list) {
    var output = []
    if (!Array.isArray(list)) return output
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object') continue
      var id = item.id === null || item.id === undefined ? '' : String(item.id)
      if (!id) continue
      output.push({
        id: id,
        name: typeof item.name === 'string' && item.name ? item.name : id,
        group: typeof item.group === 'string' ? item.group : '',
        nCells: Number(item.n_cells !== undefined ? item.n_cells : item.nCells) || 0,
      })
    }
    return output
  }

  function normalizeStringList(list) {
    if (!Array.isArray(list)) return []
    var output = []
    for (var index = 0; index < list.length; index++) {
      var text = list[index] === null || list[index] === undefined ? '' : String(list[index]).trim()
      if (text && output.indexOf(text) < 0) output.push(text)
    }
    return output
  }

  /** GET /singlecell/obs 的归一化结果。 */
  function normalizeObsPayload(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.obs && typeof payload.obs === 'object' ? payload.obs : payload)
      : {}
    var columns = []
    var list = Array.isArray(raw.columns) ? raw.columns : []
    for (var index = 0; index < list.length; index++) {
      var column = list[index]
      if (!column || typeof column !== 'object') continue
      var key = String(column.key || '')
      if (!key) continue
      columns.push({
        key: key,
        kind: column.kind === 'numeric' ? 'numeric' : 'categorical',
        nUnique: Number(column.n_unique !== undefined ? column.n_unique : column.nUnique) || 0,
        values: Array.isArray(column.values) ? column.values.map(String) : [],
      })
    }
    var groups = normalizeStringList(raw.groups)
    if (!groups.length) {
      var samplesForGroups = normalizeSamples(raw.samples)
      for (var sampleIndex = 0; sampleIndex < samplesForGroups.length; sampleIndex++) {
        var group = samplesForGroups[sampleIndex].group
        if (group && groups.indexOf(group) < 0) groups.push(group)
      }
      groups.sort()
    }
    return {
      dataset: raw.dataset === null || raw.dataset === undefined ? '' : String(raw.dataset),
      columns: columns,
      samples: normalizeSamples(raw.samples),
      groups: groups,
      groupColumn: typeof raw.group_column === 'string' && raw.group_column ? raw.group_column : 'group',
      sampleColumn: typeof raw.sample_column === 'string' && raw.sample_column ? raw.sample_column : 'sample',
    }
  }

  /** 样本表编辑行（表头：源文件名 | 细胞数 | 样本名 | 大组）。 */
  function parseSampleRows(obsPayload) {
    var obs = normalizeObsPayload(obsPayload)
    var rows = []
    for (var index = 0; index < obs.samples.length; index++) {
      var sample = obs.samples[index]
      rows.push({ id: sample.id, name: sample.name, group: sample.group, nCells: sample.nCells, dirty: false })
    }
    return rows
  }

  /** 大组下拉选项（现有 groups 优先，再补样本表里出现但列表缺的组）。 */
  function groupOptions(obsPayload) {
    var obs = normalizeObsPayload(obsPayload)
    var values = obs.groups.slice()
    for (var index = 0; index < obs.samples.length; index++) {
      var group = obs.samples[index].group
      if (group && values.indexOf(group) < 0) values.push(group)
    }
    values.sort(function (left, right) { return left < right ? -1 : (left > right ? 1 : 0) })
    var options = [{ value: '', label: '（未分组）' }]
    for (var optionIndex = 0; optionIndex < values.length; optionIndex++) {
      options.push({ value: values[optionIndex], label: values[optionIndex] })
    }
    return options
  }

  /**
   * 样本表 → POST /samples/update 的载荷。
   * 只送 id 与发生变化的 name/group 也可，但统一全量送更稳（后端按 id 覆盖）。
   */
  function buildSamplesUpdate(rows, options) {
    var settings = options && typeof options === 'object' ? options : {}
    var samples = []
    if (Array.isArray(rows)) {
      for (var index = 0; index < rows.length; index++) {
        var row = rows[index]
        if (!row || !row.id) continue
        var entry = { id: String(row.id) }
        var name = String(row.name === null || row.name === undefined ? '' : row.name).trim()
        var group = String(row.group === null || row.group === undefined ? '' : row.group).trim()
        if (name) entry.name = name
        entry.group = group
        samples.push(entry)
      }
    }
    var groups = normalizeStringList(settings.groups)
    for (var sampleIndex = 0; sampleIndex < samples.length; sampleIndex++) {
      var sampleGroup = samples[sampleIndex].group
      if (sampleGroup && groups.indexOf(sampleGroup) < 0) groups.push(sampleGroup)
    }
    var payload = { samples: samples, groups: groups }
    if (settings.sampleColumn) payload.sample_column = String(settings.sampleColumn)
    if (settings.groupColumn) payload.group_column = String(settings.groupColumn)
    return payload
  }

  /** 把 auto_group / import_map 的 mapping 合到编辑行（只改 group）。 */
  function applySampleGroupMapping(rows, mapping) {
    var map = mapping && typeof mapping === 'object' && !Array.isArray(mapping) ? mapping : {}
    var output = []
    var list = Array.isArray(rows) ? rows : []
    for (var index = 0; index < list.length; index++) {
      var row = list[index] || {}
      var next = { id: row.id, name: row.name, group: row.group, nCells: row.nCells, dirty: row.dirty === true }
      if (row.id && Object.prototype.hasOwnProperty.call(map, row.id)) {
        var group = map[row.id] === null || map[row.id] === undefined ? '' : String(map[row.id]).trim()
        if (group !== next.group) {
          next.group = group
          next.dirty = true
        }
      }
      output.push(next)
    }
    return output
  }

  function normalizeAutoGroupResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var mapping = {}
    var source = raw.mapping && typeof raw.mapping === 'object' && !Array.isArray(raw.mapping) ? raw.mapping : {}
    var keys = Object.keys(source)
    for (var index = 0; index < keys.length; index++) {
      var key = keys[index]
      var value = source[key]
      if (value === null || value === undefined) continue
      var group = String(value).trim()
      if (group) mapping[key] = group
    }
    return { mapping: mapping, unmatched: normalizeStringList(raw.unmatched) }
  }

  /** auto_group 预览表：按样本 id 排序，供"确认采用"前展示。 */
  function parseAutoGroupMapping(result) {
    var normalized = normalizeAutoGroupResult(result)
    var keys = Object.keys(normalized.mapping).sort()
    var rows = []
    for (var index = 0; index < keys.length; index++) {
      rows.push({ id: keys[index], group: normalized.mapping[keys[index]] })
    }
    return { rows: rows, unmatched: normalized.unmatched }
  }

  function normalizeImportMapResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var mapping = {}
    var list = Array.isArray(raw.updates) ? raw.updates : []
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object' || !item.id) continue
      var group = item.group === null || item.group === undefined ? '' : String(item.group).trim()
      if (group) mapping[String(item.id)] = group
    }
    var matched = Number(raw.matched) || 0
    var total = Number(raw.total) || 0
    return { updates: mapping, matched: matched, total: total }
  }

  function normalizeJoinObsResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    return {
      joinedColumns: normalizeStringList(raw.joined_columns),
      matched: Number(raw.matched) || 0,
      total: Number(raw.total) || 0,
      matchRate: Number.isFinite(Number(raw.match_rate)) ? Number(raw.match_rate) : null,
    }
  }

  /** 匹配率展示：「已匹配 850/900 细胞（94.4%）」。 */
  function formatMatchRate(result) {
    var item = result && typeof result === 'object' ? result : {}
    var matched = Number(item.matched)
    var total = Number(item.total)
    if (!Number.isFinite(matched) || !Number.isFinite(total) || total <= 0) return ''
    var rate = Number.isFinite(Number(item.matchRate))
      ? Number(item.matchRate)
      : matched / total
    return '已匹配 ' + formatNumber(matched) + '/' + formatNumber(total) + '（' + formatPct(rate) + '）'
  }

  /** 保存后的状态行：「当前分组列：group（2 组：old、young）」。 */
  function describeGroupState(groups, groupColumn) {
    var list = normalizeStringList(groups)
    var column = typeof groupColumn === 'string' && groupColumn ? groupColumn : 'group'
    if (!list.length) return '当前分组列：' + column + '（未分组）'
    return '当前分组列：' + column + '（' + list.length + ' 组：' + list.join('、') + '）'
  }

  // ==========================================================================
  // v5：数据根目录扫描 → 样本表（改名/分组/GEO 说明）→ 确认导入
  // ==========================================================================

  function basenameOf(filePath) {
    var text = String(filePath === null || filePath === undefined ? '' : filePath)
    var parts = text.split(/[\\/]/)
    return parts.length ? parts[parts.length - 1] : text
  }

  function stripDataExtension(name) {
    return String(name || '').replace(/\.(h5ad|h5|loom|rds|rdata|rda|csv|tsv|txt|mtx|gz|zip|tar|tgz|bz2)+$/i, '')
  }

  /** POST /singlecell/scan 的归一化结果。 */
  function normalizeScanResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var list = Array.isArray(raw.samples) ? raw.samples : []
    var samples = []
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object') continue
      var path = typeof item.path === 'string' ? item.path : ''
      var rawName = typeof item.raw_name === 'string' && item.raw_name
        ? item.raw_name
        : (typeof item.name === 'string' && item.name ? item.name : stripDataExtension(basenameOf(path)))
      samples.push({
        name: typeof item.name === 'string' && item.name ? item.name : rawName,
        rawName: rawName || ('样本' + (index + 1)),
        format: typeof item.format === 'string' ? item.format : '',
        formatLabel: typeof item.format_label === 'string' && item.format_label
          ? item.format_label
          : (typeof item.format === 'string' ? item.format : ''),
        path: path,
        detail: typeof item.detail === 'string' ? item.detail : '',
        trioPrefix: item.trio_prefix === null || item.trio_prefix === undefined ? '' : String(item.trio_prefix),
      })
    }
    return {
      root: typeof raw.root === 'string' ? raw.root : '',
      samples: samples,
      errors: normalizeStringList(raw.errors),
    }
  }

  /** 扫描结果 → 可编辑样本表行（默认全选、未分组）。 */
  function createScanRows(scanResult) {
    var scan = normalizeScanResult(scanResult)
    var rows = []
    var seen = {}
    for (var index = 0; index < scan.samples.length; index++) {
      var sample = scan.samples[index]
      var baseKey = (sample.path || '') + '|' + sample.trioPrefix
      var key = baseKey
      var suffix = 1
      while (seen[key]) {
        suffix += 1
        key = baseKey + '|' + suffix
      }
      seen[key] = true
      rows.push({
        key: key,
        selected: true,
        name: sample.name,
        rawName: sample.rawName,
        format: sample.format,
        formatLabel: sample.formatLabel,
        path: sample.path,
        detail: sample.detail,
        trioPrefix: sample.trioPrefix,
        group: '',
        geoTitle: '',
      })
    }
    return rows
  }

  /** 样本名重复检查（trim 后忽略大小写；返回重复名，保持首次出现写法）。 */
  function findDuplicateSampleNames(rows) {
    var counts = {}
    var display = {}
    var list = Array.isArray(rows) ? rows : []
    for (var index = 0; index < list.length; index++) {
      var name = String(list[index] && list[index].name !== undefined ? list[index].name : '').trim()
      if (!name) continue
      var key = name.toLowerCase()
      counts[key] = (counts[key] || 0) + 1
      if (!display[key]) display[key] = name
    }
    var duplicates = []
    var keys = Object.keys(counts)
    for (var keyIndex = 0; keyIndex < keys.length; keyIndex++) {
      if (counts[keys[keyIndex]] > 1) duplicates.push(display[keys[keyIndex]])
    }
    duplicates.sort()
    return duplicates
  }

  /**
   * 「确认样本与分组，开始分析」的校验：
   * 至少 1 个勾选、勾选样本都有分组、勾选样本名唯一且非空。
   */
  function validateSampleSelection(rows) {
    var list = Array.isArray(rows) ? rows : []
    var selected = []
    for (var index = 0; index < list.length; index++) {
      if (list[index] && list[index].selected === true) selected.push(list[index])
    }
    var duplicateNames = findDuplicateSampleNames(selected)
    var emptyNames = []
    var ungrouped = []
    for (var selectedIndex = 0; selectedIndex < selected.length; selectedIndex++) {
      var row = selected[selectedIndex]
      var name = String(row.name === undefined ? '' : row.name).trim()
      if (!name) emptyNames.push(row.rawName || row.path || ('第 ' + (selectedIndex + 1) + ' 个样本'))
      if (!String(row.group === undefined ? '' : row.group).trim()) {
        ungrouped.push(name || row.rawName || row.path || '')
      }
    }
    var errors = []
    if (!selected.length) errors.push('请至少勾选一个参与分析的样本')
    if (emptyNames.length) errors.push('样本名不能为空：' + emptyNames.slice(0, 5).join('、'))
    if (duplicateNames.length) errors.push('样本名重复：' + duplicateNames.join('、'))
    if (ungrouped.length) errors.push('这些样本还没有分组：' + ungrouped.slice(0, 8).join('、'))
    return {
      ok: errors.length === 0,
      errors: errors,
      selectedCount: selected.length,
      duplicateNames: duplicateNames,
      emptyNames: emptyNames,
      ungrouped: ungrouped,
    }
  }

  /** 样本表 → POST /dataset/open 的 {samples:[...]} 载荷。 */
  function buildSamplesOpenPayload(rows, options) {
    var list = Array.isArray(rows) ? rows : []
    var samples = []
    for (var index = 0; index < list.length; index++) {
      var row = list[index]
      if (!row || row.selected !== true || !row.path) continue
      var entry = {
        path: String(row.path),
        name: String(row.name === undefined ? '' : row.name).trim() || row.rawName || stripDataExtension(basenameOf(row.path)),
        group: String(row.group === undefined ? '' : row.group).trim(),
      }
      if (row.trioPrefix) entry.trio_prefix = String(row.trioPrefix)
      if (row.format) entry.format = String(row.format)
      samples.push(entry)
    }
    var payload = { samples: samples }
    var settings = options && typeof options === 'object' ? options : {}
    if (typeof settings.name === 'string' && settings.name.trim()) payload.name = settings.name.trim()
    return payload
  }

  /** 全选 / 清空 / 反选。 */
  function setScanRowsSelected(rows, selected) {
    var list = Array.isArray(rows) ? rows : []
    var value = selected === true
    return list.map(function (row) {
      return Object.assign({}, row, { selected: value })
    })
  }

  /** 批量分组：只作用于勾选行（onlySelected=false 时作用于全部）。 */
  function applyGroupToRows(rows, group, onlySelected) {
    var list = Array.isArray(rows) ? rows : []
    var value = String(group === undefined || group === null ? '' : group).trim()
    var selectedOnly = onlySelected !== false
    return list.map(function (row) {
      if (selectedOnly && row.selected !== true) return row
      return Object.assign({}, row, { group: value })
    })
  }

  /** POST /singlecell/geo 的归一化结果。 */
  function normalizeGeoResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var list = Array.isArray(raw.samples) ? raw.samples : []
    var samples = []
    for (var index = 0; index < list.length; index++) {
      var item = list[index]
      if (!item || typeof item !== 'object') continue
      var gsm = String(item.gsm || item.accession || '').trim()
      if (!gsm) continue
      samples.push({ gsm: gsm, title: String(item.title || '').trim() })
    }
    return {
      gse: String(raw.gse || '').trim(),
      title: String(raw.title || '').trim(),
      samples: samples,
    }
  }

  function sampleNameCandidates(row) {
    var values = []
    var push = function (value) {
      var text = String(value === null || value === undefined ? '' : value).trim()
      if (text && values.indexOf(text) < 0) values.push(text)
    }
    push(row.rawName)
    push(row.name)
    push(row.trioPrefix)
    if (row.path) {
      push(basenameOf(row.path))
      push(stripDataExtension(basenameOf(row.path)))
    }
    return values
  }

  /** GEO 查询结果 → 尽量把 GSM 说明匹配到扫描样本行（按文件名/trio_prefix 包含关系）。 */
  function matchGeoToSamples(rows, geoPayload) {
    var geo = normalizeGeoResult(geoPayload)
    var list = Array.isArray(rows) ? rows : []
    var matchedGsms = {}
    var outputRows = []
    for (var index = 0; index < list.length; index++) {
      var row = list[index]
      var candidates = sampleNameCandidates(row)
      var title = ''
      var gsmHit = ''
      for (var sampleIndex = 0; sampleIndex < geo.samples.length; sampleIndex++) {
        var geoSample = geo.samples[sampleIndex]
        var gsm = geoSample.gsm.toUpperCase()
        var hit = false
        for (var candidateIndex = 0; candidateIndex < candidates.length; candidateIndex++) {
          var candidate = candidates[candidateIndex].toUpperCase()
          if (candidate.indexOf(gsm) >= 0 || gsm.indexOf(candidate) >= 0) {
            hit = true
            break
          }
        }
        if (hit) {
          title = geoSample.title || geoSample.gsm
          gsmHit = geoSample.gsm
          matchedGsms[geoSample.gsm] = true
          break
        }
      }
      outputRows.push(Object.assign({}, row, { geoTitle: title }))
      if (gsmHit) gsmHit = ''
    }
    var unmatched = []
    for (var geoIndex = 0; geoIndex < geo.samples.length; geoIndex++) {
      if (!matchedGsms[geo.samples[geoIndex].gsm]) unmatched.push(geo.samples[geoIndex])
    }
    return {
      rows: outputRows,
      gse: geo.gse,
      title: geo.title,
      matched: geo.samples.length - unmatched.length,
      unmatched: unmatched,
    }
  }

  /** dataset.import_state（已跑过一半的数据）归一化。 */
  function normalizeImportState(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var steps = []
    var list = Array.isArray(raw.steps) ? raw.steps : []
    for (var index = 0; index < list.length; index++) {
      var step = String(list[index] || '')
      if (STEP_BY_ID[step] && steps.indexOf(step) < 0) steps.push(step)
    }
    var mapped = raw.mapped && typeof raw.mapped === 'object' ? raw.mapped : {}
    // 兼容已归一化过的对象（camelCase）与 Host 原始字段（snake_case）
    var pickString = function (snake, camel) {
      if (typeof mapped[snake] === 'string') return mapped[snake]
      if (typeof mapped[camel] === 'string') return mapped[camel]
      return ''
    }
    var pickNumber = function (snake, camel) {
      var value = mapped[snake] !== undefined ? mapped[snake] : mapped[camel]
      return Number(value) || 0
    }
    var pickBool = function (snake, camel) {
      return mapped[snake] === true || mapped[camel] === true
    }
    return {
      steps: steps,
      current: typeof raw.current === 'string' ? raw.current : '',
      details: raw.details && typeof raw.details === 'object' && !Array.isArray(raw.details) ? raw.details : {},
      mapped: {
        clusterColumn: pickString('cluster_column', 'clusterColumn'),
        celltypeColumn: pickString('celltype_column', 'celltypeColumn'),
        nClusters: pickNumber('n_clusters', 'nClusters'),
        nVariableFeatures: pickNumber('n_variable_features', 'nVariableFeatures'),
        hasPca: pickBool('has_pca', 'hasPca'),
        hasEmbed: pickBool('has_embed', 'hasEmbed'),
      },
    }
  }

  /** 「检测到已有进度：已完成 N 步（… 列已就绪），可继续往下跑」中文提示。 */
  function formatImportStateHint(importState, notes) {
    var state = normalizeImportState(importState)
    if (!state.steps.length && !state.mapped.celltypeColumn && !state.mapped.hasEmbed) return ''
    var parts = []
    if (state.steps.length) {
      parts.push(
        '检测到已有进度：已完成 ' + state.steps.length + ' 步（' +
        state.steps.map(stepTitleOf).join(' → ') + '），可继续往下跑',
      )
    }
    var readiness = []
    if (state.mapped.clusterColumn) {
      readiness.push('聚类列 ' + state.mapped.clusterColumn + (state.mapped.nClusters ? '（' + state.mapped.nClusters + ' 簇）' : ''))
    }
    if (state.mapped.celltypeColumn) readiness.push('细胞类型列 ' + state.mapped.celltypeColumn)
    if (state.mapped.hasEmbed) readiness.push('已有二维嵌入')
    if (state.mapped.nVariableFeatures) readiness.push('高变基因 ' + formatNumber(state.mapped.nVariableFeatures) + ' 个')
    if (readiness.length) parts.push(readiness.join('；'))
    var notesList = normalizeStringList(notes)
    var coveredPrefixes = []
    if (state.steps.length) coveredPrefixes.push('检测到已有进度')
    if (state.mapped.clusterColumn) coveredPrefixes.push('已有聚类列')
    if (state.mapped.celltypeColumn) coveredPrefixes.push('已有细胞类型列')
    if (state.mapped.hasEmbed) coveredPrefixes.push('已有二维嵌入')
    if (state.mapped.nVariableFeatures) coveredPrefixes.push('已有高变基因')
    var addedNotes = 0
    for (var index = 0; index < notesList.length && addedNotes < 2; index++) {
      var note = notesList[index]
      var covered = false
      for (var prefixIndex = 0; prefixIndex < coveredPrefixes.length; prefixIndex++) {
        if (note.indexOf(coveredPrefixes[prefixIndex]) >= 0) {
          covered = true
          break
        }
      }
      if (covered) continue
      if (parts.join('').indexOf(note) >= 0) continue
      parts.push(note)
      addedNotes += 1
    }
    return parts.join('；')
  }

  // ==========================================================================
  // v6：一键降维聚类 / 手动注释 / 分组维护 / Agent 意图解析
  // ==========================================================================

  var DIMRED_DEFAULTS = {
    engine: 'seurat',
    harmony: true,
    batchKey: '',
    resolution: 0.2,
    nPcs: 30,
    neighbors: 30,
    nTopGenes: 3000,
  }

  var DIMRED_RECIPE_ID = 'seurat6'

  function clampResolution(value) {
    var number = Number(value)
    if (!Number.isFinite(number)) return DIMRED_DEFAULTS.resolution
    return Math.round(clamp(number, 0.05, 1.5) * 100) / 100
  }

  function normalizeDimredSettings(overrides) {
    var source = overrides && typeof overrides === 'object' ? overrides : {}
    var positive = function (value, fallback) {
      var number = Number(value)
      return Number.isFinite(number) && number > 1 ? Math.floor(number) : fallback
    }
    return {
      engine: source.engine === 'builtin' ? 'builtin' : 'seurat',
      harmony: source.harmony !== false,
      batchKey: typeof source.batchKey === 'string' ? source.batchKey : '',
      resolution: clampResolution(source.resolution === undefined ? DIMRED_DEFAULTS.resolution : source.resolution),
      nPcs: positive(source.nPcs, DIMRED_DEFAULTS.nPcs),
      neighbors: positive(source.neighbors, DIMRED_DEFAULTS.neighbors),
      nTopGenes: positive(source.nTopGenes, DIMRED_DEFAULTS.nTopGenes),
    }
  }

  /**
   * 「只重聚类」的 cluster 参数。
   * Seurat 引擎：engine/use_harmony/batch_key/n_pcs/n_neighbors 都交给 Seurat 适配器；
   * 内置引擎：只发 resolution/method（保持内置 Louvain 路径）。
   */
  function buildReclusterParams(settings, autoBatchKey) {
    var normalized = normalizeDimredSettings(settings)
    var params = { resolution: normalized.resolution, method: 'louvain' }
    if (normalized.engine === 'seurat') {
      params.engine = 'seurat'
      params.use_harmony = normalized.harmony === true
      params.n_pcs = normalized.nPcs
      params.n_neighbors = normalized.neighbors
      var batchKey = normalized.batchKey || (typeof autoBatchKey === 'string' ? autoBatchKey : '')
      if (batchKey) params.batch_key = batchKey
    }
    return params
  }

  /** 「只重聚类」之后的 embed 参数（Seurat 引擎保持 Harmony 后的降维）。 */
  function buildReclusterEmbedParams(settings, autoBatchKey) {
    var normalized = normalizeDimredSettings(settings)
    var params = { layout: 'umap', n_neighbors: normalized.neighbors }
    if (normalized.engine === 'seurat') {
      params.engine = 'seurat'
      params.use_harmony = normalized.harmony === true
      params.n_pcs = normalized.nPcs
      var batchKey = normalized.batchKey || (typeof autoBatchKey === 'string' ? autoBatchKey : '')
      if (batchKey) params.batch_key = batchKey
    }
    return params
  }

  /** 「一键降维聚类」= POST /recipe {id:'seurat6'}。 */
  function buildDimredRecipeRequest(settings) {
    var normalized = normalizeDimredSettings(settings)
    return {
      id: DIMRED_RECIPE_ID,
      settings: normalized,
    }
  }

  /** 从流程卡/数据集 warnings 中挑出降维聚类相关提示（聚类数偏多/偏少等）。 */
  function pickDimredWarnings(datasetWarnings, runPayloads) {
    var output = []
    var push = function (value) {
      var text = typeof value === 'string' ? value.trim() : ''
      if (!text) return
      if (output.indexOf(text) >= 0) return
      if (output.length < 6) output.push(text)
    }
    var list = Array.isArray(datasetWarnings) ? datasetWarnings : []
    for (var index = 0; index < list.length; index++) {
      var warning = String(list[index] || '')
      if (/聚类数|resolution|harmony|Harmony|批次|主成分|近邻/i.test(warning)) push(warning)
    }
    var payloads = Array.isArray(runPayloads) ? runPayloads : []
    for (var payloadIndex = 0; payloadIndex < payloads.length; payloadIndex++) {
      var payload = payloads[payloadIndex] && typeof payloads[payloadIndex] === 'object' ? payloads[payloadIndex] : {}
      var bullets = payload.summary && Array.isArray(payload.summary.bullets) ? payload.summary.bullets : []
      for (var bulletIndex = 0; bulletIndex < bullets.length; bulletIndex++) {
        var bullet = String(bullets[bulletIndex] || '')
        if (/聚类|resolution|harmony|Harmony|近邻|主成分|cluster/i.test(bullet)) push(bullet)
      }
    }
    return output
  }

  // --------------------------------------------------------------------------
  // 手动注释工作台
  // --------------------------------------------------------------------------

  /** cluster 列表 + 已填注释 + 每簇 top 基因名（供左列展示/点击）。 */
  function buildClusterAnnotationRows(groups, annotation) {
    var list = Array.isArray(groups) ? groups : []
    var mapping = annotation && typeof annotation === 'object' ? annotation : {}
    var rows = []
    for (var index = 0; index < list.length; index++) {
      var group = list[index] || {}
      var id = group.id === null || group.id === undefined ? String(index) : String(group.id)
      var genes = Array.isArray(group.genes) ? group.genes : []
      var geneNames = []
      for (var geneIndex = 0; geneIndex < genes.length && geneNames.length < 3; geneIndex++) {
        if (genes[geneIndex] && genes[geneIndex].gene) geneNames.push(String(genes[geneIndex].gene))
      }
      rows.push({
        id: id,
        label: group.label === null || group.label === undefined ? id : String(group.label),
        nCells: Number(group.nCells) || 0,
        topGenes: geneNames,
        cellType: String(mapping[id] || '').trim(),
      })
    }
    return rows
  }

  /** 粘贴「cluster,cell_type」注释表 → {mapping, rows}。 */
  function parseAnnotationCsv(text) {
    var raw = text === null || text === undefined ? '' : String(text)
    var lines = raw.split(/\r?\n/)
    var mapping = {}
    var rows = []
    var skipped = 0
    for (var index = 0; index < lines.length; index++) {
      var line = String(lines[index] || '').trim()
      if (!line) continue
      var parts = line.split(/[,\t;，]+/).map(function (item) { return item.trim() }).filter(Boolean)
      if (parts.length < 2) {
        skipped += 1
        continue
      }
      if (index === 0 && /cluster|类型|cell.?type|label|注释/i.test(parts[1])) continue
      var cluster = parts[0]
      var label = parts.slice(1).join(' ')
      mapping[cluster] = label
      rows.push({ cluster: cluster, label: label })
    }
    return { mapping: mapping, rows: rows, skipped: skipped }
  }

  // --------------------------------------------------------------------------
  // 样本分组维护（删除空分组 / 清空分组 / 移除样本）
  // --------------------------------------------------------------------------

  function isGroupUsedBySamples(rows, group) {
    var list = Array.isArray(rows) ? rows : []
    var target = String(group === undefined || group === null ? '' : group).trim()
    if (!target) return false
    for (var index = 0; index < list.length; index++) {
      if (String(list[index] && list[index].group !== undefined ? list[index].group : '').trim() === target) return true
    }
    return false
  }

  function removeGroupFromList(groups, group) {
    var list = normalizeStringList(groups)
    var target = String(group || '')
    var output = []
    for (var index = 0; index < list.length; index++) {
      if (list[index] !== target) output.push(list[index])
    }
    return output
  }

  /** 分组下拉的特殊值「— 清空分组 —」→ 空字符串。 */
  function clearSampleGroupValue(value) {
    return value === '__clear__' ? '' : value
  }

  /** 「移除该样本」：后端暂不支持 exclude，先清空其分组。 */
  function clearSampleRowGroup(rows, sampleId) {
    var list = Array.isArray(rows) ? rows : []
    var target = String(sampleId || '')
    return list.map(function (row) {
      if (row.id !== target) return row
      return Object.assign({}, row, { group: '' })
    })
  }

  // --------------------------------------------------------------------------
  // Agent 意图解析（本地规则，纯 JS）
  // --------------------------------------------------------------------------

  function extractResolution(text) {
    var raw = String(text || '')
    var match = raw.match(/(?:resolution|分辨率|res)\s*[=:：]?\s*(\d+(?:\.\d+)?)/i)
    if (match) {
      var value = Number(match[1])
      if (Number.isFinite(value)) return clampResolution(value)
    }
    var numberMatch = raw.match(/(?:调到|改成|设为|改为|用)\s*(\d+(?:\.\d+)?)/)
    if (numberMatch) {
      var number = Number(numberMatch[1])
      if (Number.isFinite(number) && number > 0 && number <= 1.5) return clampResolution(number)
    }
    return null
  }

  var AGENT_GENE_STOPWORDS = {
    UMAP: true, PCA: true, TSV: true, CSV: true, PDF: true, HTML: true, JSON: true,
    GSE: true, GSM: true, DNA: true, RNA: true, QC: true, MT: true, RDS: true, H5AD: true,
    CDNA: true, GRN: true, CNV: true, ITH: true, GSEA: true, GSVA: true, SCENIC: true,
  }

  function extractGene(text) {
    var raw = String(text || '')
    var matches = raw.match(/\b[A-Z][A-Z0-9]{1,9}\b/g) || []
    for (var index = 0; index < matches.length; index++) {
      var token = matches[index]
      if (AGENT_GENE_STOPWORDS[token]) continue
      if (/^\d+$/.test(token)) continue
      return token
    }
    return ''
  }

  function makeAgentSuggestion(id, title, detail, action) {
    return { id: id, title: title, detail: detail, action: action }
  }

  /**
   * 用户自然语言 → 2–4 条可执行建议。
   * context: {steps:[], hasMarkers, hasAnnotation, hasGroups, sampleColumn}
   */
  function parseAgentIntent(text, context) {
    var raw = String(text === null || text === undefined ? '' : text).trim()
    var lower = raw.toLowerCase()
    var settings = context && typeof context === 'object' ? context : {}
    var steps = Array.isArray(settings.steps) ? settings.steps : []
    var matches = []
    var push = function (suggestion) {
      for (var index = 0; index < matches.length; index++) {
        if (matches[index].id === suggestion.id) return
      }
      matches.push(suggestion)
    }
    var hasResolution = extractResolution(raw)
    var gene = extractGene(raw)
    var clusterMatch = raw.match(/(?:cluster|簇)\s*[#:]?\s*(\d+)/i)
    var clusterId = clusterMatch ? String(clusterMatch[1]) : ''

    if (/harmony|整合|批次|batch/i.test(lower)) {
      push(makeAgentSuggestion(
        'harmony',
        '开启 Harmony 整合后重跑聚类',
        '用 Seurat 6.0 引擎 + Harmony 重新聚类并重算 UMAP（适合多样本/多批次数据）。',
        { type: 'recluster', params: { engine: 'seurat', use_harmony: true } },
      ))
    }
    if (hasResolution !== null) {
      push(makeAgentSuggestion(
        'resolution',
        '把 resolution 调到 ' + String(hasResolution),
        '只重跑聚类 + 嵌入，不重做 QC/归一化；参考 20–30 个 cluster。',
        { type: 'recluster', params: { resolution: hasResolution } },
      ))
    }
    if (/marker|标记基因|差异基因.*簇|top.?gene/i.test(lower)) {
      push(makeAgentSuggestion('markers', '计算各簇标记基因', '对当前 clusters 跑 t-test 标记基因（每簇 top 10）。', { type: 'markers' }))
    }
    if (/dotplot|点图/i.test(lower)) {
      push(makeAgentSuggestion('dotplot', '生成 marker DotPlot', gene ? '用基因 ' + gene + ' 生成 dotplot。' : '用各簇 top marker 生成 dotplot。', { type: 'plot', kind: 'dotplot', gene: gene }))
    }
    if (/featureplot|feature plot|表达图|画到.*umap|umap 上|着色/i.test(lower)) {
      push(makeAgentSuggestion('featureplot', '把基因表达画到 UMAP 上', gene ? '按 ' + gene + ' 表达给嵌入图着色。' : '按选中的 marker 基因给嵌入图着色。', { type: 'featureplot', gene: gene }))
    }
    if (/注释|标注|annotate|细胞类型/i.test(lower) || clusterId) {
      push(makeAgentSuggestion('annotate', '保存手动细胞类型注释', clusterId ? '把 cluster ' + clusterId + ' 的注释写入 obs.cell_type。' : '把当前注释表写入 obs.cell_type。', { type: 'annotate', cluster: clusterId }))
    }
    if (/组间差异|差异基因|de_?group|组间对比/i.test(lower)) {
      push(makeAgentSuggestion('de_group', '跑同细胞类型组间差异基因', '按分组列做逐细胞类型 Welch t + BH（需要 group 列与 clusters）。', { type: 'analysis', id: 'de_group' }))
    }
    if (/比例|proportion|占比/i.test(lower)) {
      push(makeAgentSuggestion('proportion_test', '跑组间细胞比例差异', '各组细胞类型占比 + 显著性（需要 group 列与 clusters）。', { type: 'analysis', id: 'proportion_test' }))
    }
    if (/通讯|cellchat|cell-?cell|lr\b|配体|受体/i.test(lower)) {
      push(makeAgentSuggestion('lr_communication', '跑细胞通讯分析', '内置 LR 对 + 表达乘积打分 + 置换 p（需要 clusters）。', { type: 'analysis', id: 'lr_communication' }))
    }
    if (/自检|audit|检查|qc 检查/i.test(lower)) {
      push(makeAgentSuggestion('audit', '运行自检', '对当前数据集跑确定性自检（score + 逐项问题）。', { type: 'audit' }))
    }
    if (/导出|报告|report/i.test(lower)) {
      push(makeAgentSuggestion('export_report', '导出分析报告', '导出 report（Markdown，含图与 marker 表）。', { type: 'export', what: 'report' }))
    }
    if (/methods|方法段|方法学/i.test(lower)) {
      push(makeAgentSuggestion('export_methods', '导出 Methods 段落', '导出可直接用于论文的英文 Methods（Markdown）。', { type: 'export', what: 'methods' }))
    }
    if (/一键|降维|聚类|跑流程|seurat|标准流程|继续|下一步/i.test(lower) && !/resolution|harmony/i.test(lower)) {
      push(makeAgentSuggestion('recipe_seurat6', '一键降维聚类（Seurat 6 + Harmony）', '按预设流程跑到 UMAP：HVG→PCA→Harmony→聚类→UMAP。', { type: 'recipe', id: DIMRED_RECIPE_ID }))
    }
    var keywordCount = matches.length
    if (!keywordCount) {
      return { suggestions: [], fallback: true }
    }
    // 补足 2–4 条上下文建议
    if (matches.length < 2) {
      if (steps.indexOf('cluster') < 0) {
        push(makeAgentSuggestion('recipe_seurat6', '一键降维聚类（Seurat 6 + Harmony）', '按预设流程跑到 UMAP：HVG→PCA→Harmony→聚类→UMAP。', { type: 'recipe', id: DIMRED_RECIPE_ID }))
      } else if (!settings.hasMarkers) {
        push(makeAgentSuggestion('markers', '计算各簇标记基因', '对当前 clusters 跑 t-test 标记基因（每簇 top 10）。', { type: 'markers' }))
      } else {
        push(makeAgentSuggestion('export_report', '导出分析报告', '导出 report（Markdown，含图与 marker 表）。', { type: 'export', what: 'report' }))
      }
    }
    return { suggestions: matches.slice(0, 4), fallback: false }
  }

  /** 建议动作 → 可直接发给 Agent 的中文指令。 */
  function agentActionInstruction(action) {
    var item = action && typeof action === 'object' ? action : {}
    var type = String(item.type || '')
    var params = item.params && typeof item.params === 'object' ? item.params : {}
    if (type === 'recipe') {
      return '用单细胞工作台一键降维聚类（Seurat 6.0 + Harmony，resolution ' + String(DIMRED_DEFAULTS.resolution) + '，UMAP）。'
    }
    if (type === 'recluster') {
      var pieces = []
      if (params.engine === 'seurat') pieces.push('Seurat 引擎')
      if (params.use_harmony) pieces.push('开启 Harmony')
      if (params.resolution !== undefined) pieces.push('resolution=' + String(params.resolution))
      return '用单细胞工作台只重跑聚类' + (pieces.length ? '（' + pieces.join('，') + '）' : '') + '，然后再跑二维嵌入。'
    }
    if (type === 'markers') return '用单细胞工作台计算当前聚类的标记基因（每簇 top 10）。'
    if (type === 'plot') {
      return '用单细胞工作台给当前聚类生成 ' + String(item.kind || 'dotplot') + ' 图' + (item.gene ? '（基因 ' + item.gene + '）' : '') + '。'
    }
    if (type === 'featureplot') {
      return item.gene
        ? '用单细胞工作台把基因 ' + item.gene + ' 的表达画到嵌入图上。'
        : '用单细胞工作台把选中的 marker 基因表达画到嵌入图上。'
    }
    if (type === 'annotate') {
      return '用单细胞工作台把当前手动注释表写入 obs.cell_type（cluster→细胞类型）。'
    }
    if (type === 'analysis') return '用单细胞工作台运行高级分析 ' + String(item.id || 'analysis') + '。'
    if (type === 'audit') return '用单细胞工作台对当前数据集运行自检（audit）。'
    if (type === 'export') return '用单细胞工作台导出 ' + String(item.what || 'report') + '。'
    return '用单细胞工作台继续当前数据集的分析。'
  }

  // ==========================================================================
  // API 客户端：同源相对路径 + 可注入 fetch（node 测试用假 fetch）。
  // 所有失败都转成可读中文 Error；HTTP 错误与 {ok:false,error} 都解析。
  // ==========================================================================

  function createApi(fetchImpl, baseUrl) {
    var base = typeof baseUrl === 'string' && baseUrl ? baseUrl.replace(/\/$/, '') : API

    function resolveFetch() {
      if (typeof fetchImpl === 'function') return fetchImpl
      if (typeof globalThis.fetch === 'function') return globalThis.fetch.bind(globalThis)
      return null
    }

    function apiError(message, data, status) {
      var error = new Error(message)
      error.__scApiError = true
      var payload = data && typeof data === 'object' ? data : null
      if (payload) {
        if (typeof payload.code === 'string' && payload.code) error.code = payload.code
        if (Array.isArray(payload.candidates)) error.candidates = payload.candidates
        else if (payload.data && typeof payload.data === 'object' && Array.isArray(payload.data.candidates)) {
          error.candidates = payload.data.candidates
        }
        error.data = payload
      }
      if (Number.isFinite(Number(status)) && Number(status) > 0) error.status = Number(status)
      return error
    }

    function pickError(data) {
      if (!data || typeof data !== 'object') return ''
      if (typeof data.error === 'string' && data.error) return data.error
      if (typeof data.message === 'string' && data.message) return data.message
      return ''
    }

    function parseResponse(response) {
      if (!response) return Promise.reject(apiError('Host 返回了空响应'))
      var ok = response.ok !== false
      var status = Number(response.status) || 0
      var textPromise = typeof response.text === 'function'
        ? Promise.resolve().then(function () { return response.text() }).catch(function () { return '' })
        : Promise.resolve('')
      return textPromise.then(function (text) {
        var data = null
        if (text) {
          try {
            data = JSON.parse(text)
          } catch (error) {
            data = null
          }
        }
        if (!ok) throw apiError(pickError(data) || ('HTTP ' + status), data, status)
        if (data && data.ok === false) throw apiError(pickError(data) || '请求失败', data, status)
        return data && typeof data === 'object' ? data : {}
      })
    }

    function request(path, settings, extra) {
      var doFetch = resolveFetch()
      if (!doFetch) return Promise.reject(apiError('当前环境不支持 fetch'))
      var options = Object.assign({}, settings || {})
      var extraOptions = extra && typeof extra === 'object' ? extra : {}
      var timeoutMs = Number.isFinite(extraOptions.timeoutMs) && extraOptions.timeoutMs > 0
        ? extraOptions.timeoutMs
        : DEFAULT_API_TIMEOUT_MS
      var controller = null
      var timer = null
      if (extraOptions.signal) {
        options.signal = extraOptions.signal
      } else if (typeof AbortController === 'function' && timeoutMs > 0) {
        controller = new AbortController()
        options.signal = controller.signal
        timer = setTimeout(function () {
          try {
            controller.abort()
          } catch (error) {
            /* ignore */
          }
        }, timeoutMs)
      }
      function clear() {
        if (timer) {
          clearTimeout(timer)
          timer = null
        }
      }
      var promise = Promise.resolve()
        .then(function () { return doFetch(base + path, options) })
        .then(parseResponse)
      return promise.then(function (value) {
        clear()
        return value
      }, function (error) {
        clear()
        if (error && error.__scApiError) throw error
        if (isAbortError(error)) {
          if (extraOptions.signal) throw error
          throw apiError('请求超时（' + Math.round(timeoutMs / 1000) + ' 秒）：引擎可能还在计算，请稍后刷新重试。')
        }
        throw apiError('网络请求失败：' + (formatError(error) || '无法连接 Host'))
      })
    }

    function jsonPost(path, body, extra) {
      return request(path, {
        method: 'POST',
        headers: { 'content-type': 'application/json' },
        body: JSON.stringify(body || {}),
      }, extra)
    }

    return {
      base: base,
      request: request,
      getHealth: function (extra) {
        return request('/health', { method: 'GET' }, extra)
      },
      getState: function (extra) {
        return request('/state', { method: 'GET' }, extra)
      },
      getActivity: function (extra) {
        return request('/activity', { method: 'GET' }, extra)
      },
      getFormats: function (extra) {
        return request('/formats', { method: 'GET' }, extra)
      },
      inspectSource: function (options, extra) {
        return jsonPost('/import/inspect', buildImportBody(options),
          Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      importDataset: function (options, extra) {
        return jsonPost('/dataset/open', buildImportBody(options),
          Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS * 2 }, extra || {}))
      },
      scanDataRoot: function (rootPath, recursive, extra) {
        return jsonPost('/scan', {
          path: String(rootPath || ''),
          recursive: recursive === true,
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      queryGeo: function (gse, extra) {
        return jsonPost('/geo', {
          gse: String(gse || '').trim().toUpperCase(),
          actor: 'user',
        }, Object.assign({ timeoutMs: 60000 }, extra || {}))
      },
      pickPath: function (kind, prompt, extra) {
        // Host 弹本机原生选择框（浏览器拿不到绝对路径）；timeout 留足用户操作时间
        return jsonPost('/pick', { kind: String(kind || 'dir'), prompt: String(prompt || '') },
          Object.assign({ timeoutMs: 15 * 60 * 1000 }, extra || {}))
      },
      getObs: function (datasetId, extra) {
        var query = typeof datasetId === 'string' && datasetId
          ? '?dataset=' + encodeURIComponent(datasetId)
          : ''
        return request('/obs' + query, { method: 'GET' }, extra)
      },
      samplesUpdate: function (datasetId, payload, extra) {
        var body = Object.assign({ dataset: String(datasetId || ''), actor: 'user' }, payload && typeof payload === 'object' ? payload : {})
        return jsonPost('/samples/update', body, extra)
      },
      samplesRemove: function (datasetId, samples, extra) {
        return jsonPost('/samples/remove', {
          dataset: String(datasetId || ''),
          samples: (Array.isArray(samples) ? samples : [samples]).map(function (item) { return String(item == null ? '' : item) }).filter(Boolean),
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS * 2 }, extra || {}))
      },
      samplesAutoGroup: function (datasetId, pattern, extra) {
        return jsonPost('/samples/auto-group', {
          dataset: String(datasetId || ''),
          pattern: String(pattern || ''),
          actor: 'user',
        }, extra)
      },
      samplesImportMap: function (datasetId, filePath, extra) {
        return jsonPost('/samples/import-map', {
          dataset: String(datasetId || ''),
          path: String(filePath || ''),
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      obsJoin: function (datasetId, filePath, key, extra) {
        var body = { dataset: String(datasetId || ''), path: String(filePath || ''), actor: 'user' }
        if (typeof key === 'string' && key.trim()) body.key = key.trim()
        return jsonPost('/obs/join', body, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      obsDerive: function (datasetId, options, extra) {
        var body = Object.assign({ dataset: String(datasetId || ''), actor: 'user' }, options && typeof options === 'object' ? options : {})
        return jsonPost('/obs/derive', body, extra)
      },
      obsSetGroup: function (datasetId, column, extra) {
        return jsonPost('/obs/set-group', {
          dataset: String(datasetId || ''),
          column: String(column || ''),
          actor: 'user',
        }, extra)
      },
      getAnalyses: function (datasetId, extra) {
        var query = typeof datasetId === 'string' && datasetId
          ? '?dataset=' + encodeURIComponent(datasetId)
          : ''
        return request('/analyses' + query, { method: 'GET' }, extra)
      },
      runAnalysis: function (datasetId, analysisId, params, extra) {
        return jsonPost('/analysis', {
          dataset: String(datasetId || ''),
          id: String(analysisId || ''),
          params: params && typeof params === 'object' ? params : {},
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS * 2 }, extra || {}))
      },
      openPath: function (filePath, extra) {
        return jsonPost('/dataset/open', { path: String(filePath || ''), actor: 'user' }, extra)
      },
      openDemo: function (extra) {
        return jsonPost('/dataset/demo', { actor: 'user' }, extra)
      },
      upload: function (file, filename, extra) {
        var name = typeof filename === 'string' && filename
          ? filename
          : (file && typeof file.name === 'string' && file.name ? file.name : 'upload.h5ad')
        return request('/dataset/upload', {
          method: 'POST',
          body: file,
          headers: { 'x-dsh-filename': name },
        }, extra)
      },
      selectDataset: function (datasetId, extra) {
        return jsonPost('/dataset/select', { dataset: String(datasetId || '') }, extra)
      },
      closeDataset: function (datasetId, extra) {
        return jsonPost('/dataset/close', { dataset: String(datasetId || ''), actor: 'user' }, extra)
      },
      preview: function (datasetId, extra) {
        return jsonPost('/preview', { dataset: String(datasetId || '') }, extra)
      },
      runStep: function (datasetId, step, params, extra) {
        return jsonPost('/run', {
          dataset: String(datasetId || ''),
          step: String(step || ''),
          params: params && typeof params === 'object' ? params : {},
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      getEmbedding: function (datasetId, colorBy, extra) {
        return jsonPost('/embedding', {
          dataset: String(datasetId || ''),
          color_by: colorBy && typeof colorBy === 'object' ? colorBy : { kind: 'cluster', key: 'clusters' },
        }, extra)
      },
      plot: function (datasetId, kind, params, extra) {
        return jsonPost('/plot', {
          dataset: String(datasetId || ''),
          kind: String(kind || ''),
          params: params && typeof params === 'object' ? params : {},
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      markers: function (datasetId, options, extra) {
        var settings = options && typeof options === 'object' ? options : {}
        return jsonPost('/markers', {
          dataset: String(datasetId || ''),
          top_n: Number(settings.top_n) > 0 ? Number(settings.top_n) : 10,
          method: typeof settings.method === 'string' && settings.method ? settings.method : 't-test',
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      getAudit: function (datasetId, extra) {
        return jsonPost('/audit', { dataset: String(datasetId || ''), actor: 'user' }, extra)
      },
      getAutoannotate: function (datasetId, extra) {
        return jsonPost('/autoannotate', { dataset: String(datasetId || ''), actor: 'user' },
          Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      runRecipe: function (datasetId, recipeId, extra) {
        return jsonPost('/recipe', {
          dataset: String(datasetId || ''),
          id: String(recipeId || ''),
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS * 2 }, extra || {}))
      },
      resolutionScan: function (datasetId, resolutions, extra) {
        var list = Array.isArray(resolutions)
          ? resolutions
          : parseResolutions(resolutions).list
        return jsonPost('/resolution-scan', {
          dataset: String(datasetId || ''),
          resolutions: list,
          actor: 'user',
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      exportData: function (datasetId, what, path, extra) {
        var body = { dataset: String(datasetId || ''), what: String(what || 'h5ad'), actor: 'user' }
        if (typeof path === 'string' && path.trim()) body.path = path.trim()
        return jsonPost('/export', body, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      reset: function (extra) {
        return jsonPost('/reset', { actor: 'user' }, extra)
      },
    }
  }

  // ==========================================================================
  // i18n：用宿主 locale 服务造 t()；取不到时回落中文词典，保证不出现裸 key。
  // ==========================================================================

  function createTranslator(locale) {
    var bound = null
    try {
      if (locale && typeof locale.bind === 'function') bound = locale.bind(NS)
    } catch (error) {
      bound = null
    }
    return function translate(key) {
      if (typeof bound === 'function') {
        try {
          var value = bound(key)
          if (typeof value === 'string' && value && value !== key) return value
        } catch (error) {
          /* 宿主 locale 不可用时静默回落 */
        }
      }
      var fallback = DICT_ZH[key]
      return typeof fallback === 'string' ? fallback : String(key)
    }
  }

  /** 模板取词：字典里用 `{v}` 占位，避免依赖宿主 t() 的插值实现。 */
  function tf(translate, key, value) {
    var raw = typeof translate === 'function' ? translate(key) : ''
    var template = typeof raw === 'string' && raw && raw !== key ? raw : (DICT_ZH[key] || '')
    return String(template).replace(/\{v\}/g, value === null || value === undefined ? '' : String(value))
  }

  // ==========================================================================
  // 纯函数测试出口。node 环境（typeof window === 'undefined'）挂到
  // globalThis.__SC_TEST__；浏览器里也会挂进 factory 返回的 __test__，便于排查。
  // ==========================================================================
  var SC_TEST = {
    NS: NS,
    API: API,
    POLL_MS: POLL_MS,
    MAX_HISTORY_ROWS: MAX_HISTORY_ROWS,
    PALETTE: PALETTE,
    CONTINUOUS_STOPS: CONTINUOUS_STOPS,
    AUDIT_LEVEL_STYLE: AUDIT_LEVEL_STYLE,
    EXPORT_KINDS: EXPORT_KINDS,
    STEPS: STEPS,
    STEP_BY_ID: STEP_BY_ID,
    STEP_TITLES_ZH: STEP_TITLES_ZH,
    DICT_ZH: DICT_ZH,
    DICT_EN: DICT_EN,
    clamp: clamp,
    formatNumber: formatNumber,
    formatPct: formatPct,
    formatPValue: formatPValue,
    formatDuration: formatDuration,
    formatClock: formatClock,
    formatError: formatError,
    plainText: plainText,
    isAbortError: isAbortError,
    stepTitleOf: stepTitleOf,
    doneStepSet: doneStepSet,
    isStepReady: isStepReady,
    nextRunnableSteps: nextRunnableSteps,
    nextStep: nextStep,
    stepStatusOf: stepStatusOf,
    stepStatusKey: stepStatusKey,
    defaultParams: defaultParams,
    collectStepParams: collectStepParams,
    formatParamValue: formatParamValue,
    formatParamList: formatParamList,
    buildAgentInstruction: buildAgentInstruction,
    createView: createView,
    viewCenter: viewCenter,
    worldToScreen: worldToScreen,
    screenToWorld: screenToWorld,
    zoomViewAt: zoomViewAt,
    panView: panView,
    computeBounds: computeBounds,
    fitView: fitView,
    findNearestPoint: findNearestPoint,
    tooltipText: tooltipText,
    hexToRgb: hexToRgb,
    lerpColor: lerpColor,
    pickContinuousColor: pickContinuousColor,
    normalizeEngine: normalizeEngine,
    chooseEngine: chooseEngine,
    engineBadge: engineBadge,
    normalizeStepIds: normalizeStepIds,
    normalizeHistory: normalizeHistory,
    normalizeDataset: normalizeDataset,
    normalizeSuggestions: normalizeSuggestions,
    normalizeActivity: normalizeActivity,
    normalizeRecipes: normalizeRecipes,
    normalizeRecipeProgress: normalizeRecipeProgress,
    normalizeState: normalizeState,
    findDataset: findDataset,
    normalizePreview: normalizePreview,
    normalizeCategory: normalizeCategory,
    normalizeEmbedding: normalizeEmbedding,
    normalizeMarkersGroups: normalizeMarkersGroups,
    normalizeFigure: normalizeFigure,
    normalizeAudit: normalizeAudit,
    normalizeAutoannotate: normalizeAutoannotate,
    normalizeResolutionScan: normalizeResolutionScan,
    normalizeExportResult: normalizeExportResult,
    formatHistoryEntry: formatHistoryEntry,
    buildAnnotationMapping: buildAnnotationMapping,
    parseResolutions: parseResolutions,
    recipeProgressText: recipeProgressText,
    activityIcon: activityIcon,
    createApi: createApi,
    createTranslator: createTranslator,
    DIMRED_DEFAULTS: DIMRED_DEFAULTS,
    DIMRED_RECIPE_ID: DIMRED_RECIPE_ID,
    clampResolution: clampResolution,
    normalizeDimredSettings: normalizeDimredSettings,
    buildReclusterParams: buildReclusterParams,
    buildReclusterEmbedParams: buildReclusterEmbedParams,
    buildDimredRecipeRequest: buildDimredRecipeRequest,
    pickDimredWarnings: pickDimredWarnings,
    buildClusterAnnotationRows: buildClusterAnnotationRows,
    parseAnnotationCsv: parseAnnotationCsv,
    isGroupUsedBySamples: isGroupUsedBySamples,
    removeGroupFromList: removeGroupFromList,
    clearSampleGroupValue: clearSampleGroupValue,
    clearSampleRowGroup: clearSampleRowGroup,
    extractResolution: extractResolution,
    extractGene: extractGene,
    parseAgentIntent: parseAgentIntent,
    agentActionInstruction: agentActionInstruction,
    basenameOf: basenameOf,
    stripDataExtension: stripDataExtension,
    normalizeScanResult: normalizeScanResult,
    createScanRows: createScanRows,
    findDuplicateSampleNames: findDuplicateSampleNames,
    validateSampleSelection: validateSampleSelection,
    buildSamplesOpenPayload: buildSamplesOpenPayload,
    setScanRowsSelected: setScanRowsSelected,
    applyGroupToRows: applyGroupToRows,
    normalizeGeoResult: normalizeGeoResult,
    matchGeoToSamples: matchGeoToSamples,
    normalizeImportState: normalizeImportState,
    formatImportStateHint: formatImportStateHint,
    parsePathLines: parsePathLines,
    normalizeSamples: normalizeSamples,
    normalizeObsPayload: normalizeObsPayload,
    parseSampleRows: parseSampleRows,
    groupOptions: groupOptions,
    buildSamplesUpdate: buildSamplesUpdate,
    applySampleGroupMapping: applySampleGroupMapping,
    normalizeAutoGroupResult: normalizeAutoGroupResult,
    parseAutoGroupMapping: parseAutoGroupMapping,
    normalizeImportMapResult: normalizeImportMapResult,
    normalizeJoinObsResult: normalizeJoinObsResult,
    formatMatchRate: formatMatchRate,
    describeGroupState: describeGroupState,
    formatBytes: formatBytes,
    normalizeFormats: normalizeFormats,
    normalizeObsColumns: normalizeObsColumns,
    normalizeImportInspect: normalizeImportInspect,
    normalizeCandidate: normalizeCandidate,
    formatCandidates: formatCandidates,
    buildImportBody: buildImportBody,
    ANALYSIS_GROUPS: ANALYSIS_GROUPS,
    ANALYSIS_GROUP_LABEL_KEYS: ANALYSIS_GROUP_LABEL_KEYS,
    ANALYSIS_GROUP_OTHER: ANALYSIS_GROUP_OTHER,
    IMPORT_ACCEPT: IMPORT_ACCEPT,
    FALLBACK_FORMATS: FALLBACK_FORMATS,
    normalizeAnalysisParam: normalizeAnalysisParam,
    normalizeAnalysisParams: normalizeAnalysisParams,
    normalizeAnalyses: normalizeAnalyses,
    normalizeAnalysesPayload: normalizeAnalysesPayload,
    mergeObsColumns: mergeObsColumns,
    groupAnalyses: groupAnalyses,
    paramOptionValues: paramOptionValues,
    paramDefaultValue: paramDefaultValue,
    paramDefaults: paramDefaults,
    collectAnalysisParams: collectAnalysisParams,
    needsHint: needsHint,
    formatTableCell: formatTableCell,
    formatTableRows: formatTableRows,
    normalizeAnalysisTable: normalizeAnalysisTable,
    normalizeAnalysisFigure: normalizeAnalysisFigure,
    normalizeAnalysisResult: normalizeAnalysisResult,
    FIXED_COLOR_OPTIONS: FIXED_COLOR_OPTIONS,
    buildColorOptions: buildColorOptions,
    tf: tf,
  }

  // ==========================================================================
  // 浏览器层：整个 factory 只在宿主物化模块时执行；node 下永远不会被调用。
  // ==========================================================================

  function factory(require) {
    var React = require('react')

    /** React.createElement 简写（不写 JSX，避免宿主构建链依赖）。 */
    function h(type, props) {
      var children = Array.prototype.slice.call(arguments, 2)
      return React.createElement.apply(React, [type, props].concat(children))
    }

    // t()：apply() 里用宿主 locale 服务重新绑定；默认回落中文词典。
    var T = function (key) {
      var fallback = DICT_ZH[key]
      return typeof fallback === 'string' ? fallback : String(key)
    }

    function warnPlugin(message, error) {
      try {
        console.warn('[dsh-singlecell] ' + message, error || '')
      } catch (ignored) {
        /* console 不可用就算了 */
      }
    }

    function summaryTextOf(summary) {
      if (!summary || typeof summary !== 'object') return ''
      if (Array.isArray(summary.bullets) && summary.bullets.length) return plainText(summary.bullets.join('；'), 220)
      if (typeof summary.title === 'string') return plainText(summary.title, 220)
      return ''
    }

    // ==========================================================================
    // 样式：一个 CSS 字符串，由 <style> 渲染在组件内部（生命周期跟随组件）。
    // 颜色只用 --dsw-alias-* 主题变量，派生色统一走 color-mix，深/浅色自适应。
    // 所有类名 sc- 前缀，避免和宿主/其它插件撞车。
    // ==========================================================================
    var CSS = [
      '.sc-root{display:flex;flex-direction:column;height:100%;min-height:72vh;overflow:hidden;position:relative;isolation:isolate;color:var(--dsw-alias-label-primary);background:var(--dsw-alias-bg-base);border-radius:12px;font-size:13px;line-height:1.5}',
      '.sc-root *,.sc-root *::before,.sc-root *::after{box-sizing:border-box}',
      '.sc-root button,.sc-root input,.sc-root select,.sc-root textarea{font-family:inherit}',

      // 顶栏
      '.sc-top{display:flex;align-items:center;gap:10px;flex-wrap:wrap;flex:none;min-height:50px;padding:8px 12px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sc-top-main{display:flex;align-items:center;gap:10px;flex-wrap:wrap;min-width:0;flex:1 1 320px}',
      '.sc-title-block{display:flex;flex-direction:column;gap:0;min-width:0}',
      '.sc-title{font-size:14px;font-weight:600;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sc-subtitle{font-size:11px;color:var(--dsw-alias-label-secondary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sc-badge{display:inline-flex;align-items:center;gap:5px;height:21px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-size:11px;line-height:1;cursor:default;white-space:nowrap}',
      '.sc-badge-dot{width:6px;height:6px;border-radius:50%;flex:none;background:currentColor;opacity:.9}',
      '.sc-badge-ok{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 10%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-success-primary)}',
      '.sc-badge-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 40%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 10%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-warn-primary)}',
      '.sc-badge-bad{border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 10%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-error-primary)}',
      '.sc-chips{display:flex;align-items:center;gap:5px;flex-wrap:wrap;min-width:0}',
      '.sc-chip{display:inline-flex;align-items:center;gap:4px;max-width:180px;height:23px;padding:0 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary);font:inherit;font-size:11.5px;cursor:pointer;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sc-chip:hover:not(:disabled){border-color:var(--dsw-alias-border-l2)}',
      '.sc-chip-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 10%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',
      '.sc-top-actions{display:flex;align-items:center;gap:6px;flex-wrap:wrap;margin-left:auto;justify-content:flex-end}',
      '.sc-open-group{display:flex;align-items:center;gap:5px;min-width:0}',
      '.sc-form-label{font-size:11px;color:var(--dsw-alias-label-tertiary);white-space:nowrap}',
      '.sc-path-input{width:230px;max-width:32vw}',
      '.sc-recipe{display:flex;align-items:center;gap:5px}',
      '.sc-recipe-select{max-width:170px}',

      // 表单控件
      '.sc-btn{display:inline-flex;align-items:center;justify-content:center;gap:4px;height:28px;padding:0 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary);font:inherit;font-size:12px;line-height:1;cursor:pointer;white-space:nowrap;transition:border-color .12s ease,background .12s ease,opacity .12s ease}',
      '.sc-btn:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-overlay)}',
      '.sc-btn:disabled{opacity:.45;cursor:default}',
      '.sc-btn-primary{background:var(--dsw-alias-brand-primary);border-color:var(--dsw-alias-brand-primary);color:var(--dsw-alias-bg-base)}',
      '.sc-btn-primary:hover:not(:disabled){opacity:.88;background:var(--dsw-alias-brand-primary);border-color:var(--dsw-alias-brand-primary)}',
      '.sc-btn-sm{height:24px;padding:0 8px;font-size:11.5px;border-radius:6px}',
      '.sc-input{height:28px;padding:0 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);color:var(--dsw-alias-label-primary);font:inherit;font-size:12px;outline:none;min-width:0}',
      '.sc-input:focus{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1))}',
      '.sc-input:disabled{opacity:.55}',
      '.sc-select{cursor:pointer;padding-right:4px}',

      // 三栏
      '.sc-body{display:flex;flex:1 1 auto;min-height:0;overflow:hidden}',
      '.sc-left{display:flex;flex-direction:column;flex:0 0 268px;width:268px;min-height:0;border-right:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sc-right{display:flex;flex-direction:column;flex:0 0 300px;width:300px;min-height:0;overflow-y:auto;border-left:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sc-center{display:flex;flex-direction:column;flex:1 1 auto;min-width:0;min-height:0;background:var(--dsw-alias-bg-base)}',
      '.sc-panel-head{display:flex;align-items:center;justify-content:space-between;gap:8px;flex:none;padding:8px 12px;border-bottom:1px solid var(--dsw-alias-border-l1);font-size:11px;font-weight:600;letter-spacing:.03em;color:var(--dsw-alias-label-secondary);text-transform:none}',
      '.sc-panel-title{color:var(--dsw-alias-label-secondary)}',
      '.sc-panel-note{color:var(--dsw-alias-label-tertiary);font-weight:400}',
      '.sc-panel-hint{padding:8px 12px 0;font-size:11px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-muted{font-size:11.5px;color:var(--dsw-alias-label-secondary)}',
      '.sc-section-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11px;font-weight:600;letter-spacing:.03em;color:var(--dsw-alias-label-secondary);margin-top:2px}',
      '.sc-inline-err{padding:5px 8px;border-radius:6px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,transparent);background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 8%,transparent);color:var(--dsw-alias-state-error-primary);font-size:11.5px;word-break:break-word}',
      '.sc-ok{color:var(--dsw-alias-state-success-primary);font-size:11.5px}',
      '.sc-bad{color:var(--dsw-alias-state-error-primary);font-size:11px}',
      '.sc-warn-text{color:var(--dsw-alias-state-warn-primary);font-size:11.5px}',
      '.sc-divider{flex:none;height:1px;margin:2px 0;background:var(--dsw-alias-border-l1)}',
      '.sc-inline-row{display:flex;align-items:center;gap:6px;flex-wrap:wrap;min-width:0}',

      // 左栏：流程卡片
      '.sc-steps{display:flex;flex-direction:column;gap:8px;padding:10px 12px 16px;overflow-y:auto;min-height:0}',
      '.sc-step{display:flex;flex-direction:column;gap:5px;padding:8px 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-base)}',
      '.sc-step-done{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 28%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 5%,var(--dsw-alias-bg-base))}',
      '.sc-step-running{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,var(--dsw-alias-bg-base))}',
      '.sc-step-blocked{opacity:.78}',
      '.sc-step-head{display:flex;align-items:center;gap:7px;cursor:pointer;min-width:0}',
      '.sc-step-head:focus-visible{outline:2px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 55%,transparent);outline-offset:2px;border-radius:6px}',
      '.sc-step-dot{display:inline-flex;align-items:center;justify-content:center;flex:none;width:20px;height:20px;border-radius:50%;border:1px solid var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-size:10.5px;font-weight:600}',
      '.sc-step-dot-done{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 55%,transparent);background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 14%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-success-primary)}',
      '.sc-step-dot-running{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 55%,transparent);background:color-mix(in srgb,var(--dsw-alias-brand-primary) 14%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',
      '.sc-step-dot-pending{border-color:var(--dsw-alias-border-l2);color:var(--dsw-alias-label-secondary)}',
      '.sc-step-dot-blocked{border-color:var(--dsw-alias-border-l1);color:var(--dsw-alias-label-tertiary)}',
      '.sc-step-title{flex:1 1 auto;min-width:0;font-size:12.5px;font-weight:500;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sc-step-status{flex:none;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-status-done{color:var(--dsw-alias-state-success-primary)}',
      '.sc-status-running{color:var(--dsw-alias-brand-primary)}',
      '.sc-status-blocked{color:var(--dsw-alias-label-tertiary)}',
      '.sc-step-caret{flex:none;font-size:10px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-step-params{display:flex;flex-direction:column;gap:6px;padding:6px 2px 2px;border-top:1px dashed var(--dsw-alias-border-l1);margin-top:2px}',
      '.sc-param{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11.5px;color:var(--dsw-alias-label-secondary)}',
      '.sc-param-label{flex:0 1 auto;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-param .sc-input{flex:0 0 108px;width:108px;height:24px;font-size:11.5px;padding:0 6px}',
      '.sc-param-hint{font-size:11px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-step-foot{display:flex;align-items:center;gap:6px;flex-wrap:wrap}',
      '.sc-step-summary{flex:1 1 auto;min-width:0;font-size:11px;color:var(--dsw-alias-label-tertiary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sc-step .sc-inline-err{margin-top:2px}',

      // 中栏 tabs
      '.sc-tabs{display:flex;align-items:center;gap:2px;flex:none;padding:6px 10px 0;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sc-tab{display:inline-flex;align-items:center;justify-content:center;gap:5px;height:30px;padding:0 12px;border:none;border-bottom:2px solid transparent;border-radius:7px 7px 0 0;background:transparent;color:var(--dsw-alias-label-secondary);font:inherit;font-size:12px;cursor:pointer;white-space:nowrap}',
      '.sc-tab:hover{color:var(--dsw-alias-label-primary);background:var(--dsw-alias-bg-layer-2)}',
      '.sc-tab-on{color:var(--dsw-alias-brand-primary);border-bottom-color:var(--dsw-alias-brand-primary);font-weight:600}',
      '.sc-pane{display:flex;flex:1 1 auto;flex-direction:column;gap:10px;min-height:0;padding:10px 12px 16px;overflow-y:auto}',
      '.sc-toolbar{display:flex;align-items:center;gap:8px;flex-wrap:wrap;flex:none;min-height:28px}',
      '.sc-canvas-stack{display:flex;flex-direction:column;flex:1 1 auto;gap:8px;min-height:300px;min-width:0}',
      '.sc-canvas-wrap{position:relative;flex:1 1 auto;min-height:300px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-base);overflow:hidden}',
      '.sc-canvas{display:block;width:100%;height:100%;cursor:grab;touch-action:none}',
      '.sc-canvas:active{cursor:grabbing}',
      '.sc-tip{position:absolute;z-index:3;pointer-events:none;max-width:220px;padding:4px 8px;border:1px solid var(--dsw-alias-border-l2);border-radius:6px;background:var(--dsw-alias-bg-overlay);color:var(--dsw-alias-label-primary);font-size:11px;white-space:nowrap;box-shadow:0 4px 14px color-mix(in srgb,var(--dsw-alias-label-primary) 16%,transparent)}',
      '.sc-canvas-meta{position:absolute;left:8px;bottom:6px;z-index:2;font-size:10.5px;color:var(--dsw-alias-label-tertiary);pointer-events:none}',
      '.sc-canvas-reset{position:absolute;right:8px;top:8px;z-index:2}',
      '.sc-embed-tools{display:flex;align-items:center;gap:8px;flex-wrap:wrap}',
      '.sc-tool-field{display:inline-flex;align-items:center;gap:6px;font-size:11.5px;color:var(--dsw-alias-label-secondary)}',
      '.sc-gene-form{display:inline-flex;align-items:center;gap:5px}',
      '.sc-gene-input{width:150px}',
      '.sc-legend{display:flex;flex-direction:column;gap:5px;flex:none;max-height:110px;overflow-y:auto}',
      '.sc-legend-row{display:flex;flex-wrap:wrap;gap:5px}',
      '.sc-legend-chip{display:inline-flex;align-items:center;gap:5px;height:22px;padding:0 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary);font:inherit;font-size:11px;cursor:pointer;white-space:nowrap}',
      '.sc-legend-chip:hover{border-color:var(--dsw-alias-border-l2)}',
      '.sc-legend-chip-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 55%,var(--dsw-alias-border-l1));box-shadow:0 0 0 1px color-mix(in srgb,var(--dsw-alias-brand-primary) 35%,transparent)}',
      '.sc-legend-dot{width:8px;height:8px;border-radius:50%;flex:none}',
      '.sc-legend-n{color:var(--dsw-alias-label-tertiary);font-size:10.5px}',
      '.sc-gradient{display:flex;align-items:center;gap:8px;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-gradient-bar{flex:1 1 auto;height:10px;max-width:280px;border-radius:999px;border:1px solid var(--dsw-alias-border-l1)}',
      '.sc-notice{display:flex;align-items:center;gap:6px;padding:10px 12px;border:1px dashed var(--dsw-alias-border-l2);border-radius:10px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-secondary);font-size:12px}',
      '.sc-guide{display:flex;flex-direction:column;align-items:flex-start;gap:8px;flex:1 1 auto;justify-content:center;padding:22px 20px;border:1.5px dashed var(--dsw-alias-border-l2);border-radius:12px;background:var(--dsw-alias-bg-layer-1);text-align:left}',
      '.sc-guide-title{font-size:13.5px;font-weight:600}',
      '.sc-guide-hint{font-size:12px;color:var(--dsw-alias-label-secondary);max-width:460px}',
      '.sc-hint-line{font-size:11.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-check{display:inline-flex;align-items:center;gap:5px;font-size:12px;color:var(--dsw-alias-label-secondary);cursor:pointer;user-select:none}',
      '.sc-check input{accent-color:var(--dsw-alias-brand-primary)}',
      '.sc-figure{display:flex;flex-direction:column;gap:4px;align-items:flex-start;max-width:100%}',
      '.sc-figure img{display:block;max-width:100%;height:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-2)}',
      '.sc-figure-meta{font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-figure-sm img{max-width:520px}',
      '.sc-table-wrap{flex:1 1 auto;min-height:120px;overflow:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.sc-table{width:100%;border-collapse:collapse;font-size:11.5px}',
      '.sc-table th{position:sticky;top:0;z-index:1;padding:6px 10px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-secondary);font-weight:600;text-align:left;white-space:nowrap}',
      '.sc-table td{padding:4px 10px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 55%,transparent);color:var(--dsw-alias-label-primary);white-space:nowrap}',
      '.sc-table tbody tr:last-child td{border-bottom:none}',
      '.sc-td-group{vertical-align:top;font-weight:600;color:var(--dsw-alias-brand-primary);background:color-mix(in srgb,var(--dsw-alias-brand-primary) 5%,transparent)}',
      '.sc-td-num{text-align:right;font-variant-numeric:tabular-nums}',

      // 聚类 tab
      '.sc-cluster-list{display:flex;flex-direction:column;gap:5px;max-height:320px;overflow-y:auto;padding-right:2px}',
      '.sc-cluster-row{display:flex;align-items:center;gap:8px;padding:4px 6px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-1)}',
      '.sc-cluster-id{flex:0 0 auto;min-width:34px;font-weight:600;font-size:12px;color:var(--dsw-alias-brand-primary)}',
      '.sc-cluster-n{flex:0 0 auto;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-cluster-input{flex:1 1 auto;min-width:0;height:24px;font-size:11.5px}',

      // 日志 tab
      '.sc-log-list{display:flex;flex-direction:column;gap:4px}',
      '.sc-log-row{display:grid;grid-template-columns:66px 20px minmax(90px,150px) 1fr auto;align-items:baseline;gap:8px;padding:4px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-1);font-size:11.5px}',
      '.sc-log-time{color:var(--dsw-alias-label-tertiary);font-variant-numeric:tabular-nums;font-size:10.5px}',
      '.sc-log-actor{font-size:12px;text-align:center}',
      '.sc-log-action{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:11px;color:var(--dsw-alias-brand-primary);overflow:hidden;text-overflow:ellipsis}',
      '.sc-log-detail{color:var(--dsw-alias-label-secondary);overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-log-ms{color:var(--dsw-alias-label-tertiary);font-size:10.5px;white-space:nowrap}',

      // 右栏 Agent
      '.sc-agent{display:flex;flex-direction:column;gap:8px;padding:10px 12px 14px;min-height:0}',
      '.sc-agent-status{display:flex;flex-direction:column;gap:3px;padding:9px 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-base)}',
      '.sc-agent-name{font-size:13px;font-weight:600;word-break:break-word}',
      '.sc-agent-meta{font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-agent-row{display:flex;align-items:center;justify-content:space-between;gap:8px}',
      '.sc-agent-k{font-size:11px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-suggest{display:flex;flex-direction:column;gap:5px}',
      '.sc-sug{display:flex;flex-direction:column;gap:4px;padding:8px 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-base)}',
      '.sc-sug-title{font-size:12px;font-weight:500}',
      '.sc-sug-detail{font-size:11px;color:var(--dsw-alias-label-secondary);word-break:break-word}',
      '.sc-sug-actions{display:flex;justify-content:flex-end}',
      '.sc-timeline{display:flex;flex-direction:column;gap:2px}',
      '.sc-act{display:grid;grid-template-columns:16px 1fr;gap:5px;padding:3px 2px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 50%,transparent)}',
      '.sc-act:last-child{border-bottom:none}',
      '.sc-act-icon{font-size:11px;line-height:1.6;text-align:center}',
      '.sc-act-line{display:flex;align-items:baseline;gap:6px;flex-wrap:wrap;min-width:0}',
      '.sc-act-action{font-size:11px;color:var(--dsw-alias-label-primary);overflow:hidden;text-overflow:ellipsis}',
      '.sc-act-detail{font-size:10.5px;color:var(--dsw-alias-label-secondary);word-break:break-word}',
      '.sc-act-meta{font-size:10px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-agent-hint{margin-top:auto;padding:8px 10px;border:1px dashed var(--dsw-alias-border-l1);border-radius:8px;font-size:11px;color:var(--dsw-alias-label-tertiary)}',

      // 记录 tab（v2）
      '.sc-record{display:flex;flex-direction:column;gap:12px}',
      '.sc-record-block{display:flex;flex-direction:column;gap:6px;padding:9px 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.sc-record-head{display:flex;align-items:center;justify-content:space-between;gap:8px;flex-wrap:wrap}',
      '.sc-record-title{font-size:12px;font-weight:600}',
      '.sc-record-note{font-size:11px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-history{display:flex;flex-direction:column;gap:3px;max-height:260px;overflow-y:auto}',
      '.sc-history-row{display:grid;grid-template-columns:62px minmax(70px,110px) minmax(110px,1fr) minmax(120px,1.4fr) auto;gap:8px;align-items:baseline;padding:4px 6px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 55%,transparent);font-size:11.5px}',
      '.sc-history-row:last-child{border-bottom:none}',
      '.sc-history-head{color:var(--dsw-alias-label-tertiary);font-size:10.5px;border-bottom:1px solid var(--dsw-alias-border-l1)}',
      '.sc-history-time{color:var(--dsw-alias-label-tertiary);font-size:10.5px;font-variant-numeric:tabular-nums}',
      '.sc-history-step{font-weight:500;color:var(--dsw-alias-label-primary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sc-history-params{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:10.5px;color:var(--dsw-alias-label-secondary);overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-history-summary{color:var(--dsw-alias-label-secondary);overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-history-ms{color:var(--dsw-alias-label-tertiary);font-size:10.5px;white-space:nowrap;text-align:right}',
      '.sc-audit-head{display:flex;align-items:center;gap:8px;flex-wrap:wrap}',
      '.sc-audit-score{display:inline-flex;align-items:baseline;gap:5px;padding:2px 10px;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-2);font-size:12px;font-weight:600}',
      '.sc-audit-score-ok{color:var(--dsw-alias-state-success-primary);border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1))}',
      '.sc-audit-score-warn{color:var(--dsw-alias-state-warn-primary);border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 40%,var(--dsw-alias-border-l1))}',
      '.sc-audit-score-fail{color:var(--dsw-alias-state-error-primary);border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,var(--dsw-alias-border-l1))}',
      '.sc-check-list{display:flex;flex-direction:column;gap:4px}',
      '.sc-check-row{display:grid;grid-template-columns:18px 1fr;gap:6px;padding:5px 7px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base)}',
      '.sc-check-row-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 40%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 7%,var(--dsw-alias-bg-base))}',
      '.sc-check-row-fail{border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 7%,var(--dsw-alias-bg-base))}',
      '.sc-check-icon{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;border-radius:50%;font-size:10px;font-weight:700;border:1px solid var(--dsw-alias-border-l2);color:var(--dsw-alias-label-secondary)}',
      '.sc-check-icon-ok{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 55%,transparent);color:var(--dsw-alias-state-success-primary)}',
      '.sc-check-icon-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 55%,transparent);color:var(--dsw-alias-state-warn-primary)}',
      '.sc-check-icon-fail{border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 55%,transparent);color:var(--dsw-alias-state-error-primary)}',
      '.sc-check-title{font-size:11.5px;font-weight:500}',
      '.sc-check-detail{font-size:11px;color:var(--dsw-alias-label-secondary);word-break:break-word}',
      '.sc-check-hint{font-size:10.5px;color:var(--dsw-alias-label-tertiary);word-break:break-word}',
      '.sc-ann-list{display:flex;flex-direction:column;gap:4px;max-height:300px;overflow-y:auto}',
      '.sc-ann-row{display:grid;grid-template-columns:34px minmax(70px,110px) 1fr auto;gap:8px;align-items:center;padding:5px 7px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base)}',
      '.sc-ann-cluster{font-weight:600;font-size:11.5px;color:var(--dsw-alias-brand-primary)}',
      '.sc-ann-label{font-size:11.5px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-ann-meta{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.sc-ann-evidence{font-size:10.5px;color:var(--dsw-alias-label-tertiary);overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-ann-actions{display:flex;gap:4px;flex-wrap:wrap;justify-content:flex-end}',
      '.sc-scan-form{display:flex;align-items:center;gap:6px;flex-wrap:wrap}',
      '.sc-scan-input{width:180px}',
      '.sc-export-btns{display:flex;align-items:center;gap:5px;flex-wrap:wrap}',
      '.sc-export-result{display:flex;flex-direction:column;gap:3px;padding:6px 8px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-success-primary) 35%,var(--dsw-alias-border-l1));border-radius:8px;background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 6%,var(--dsw-alias-bg-base))}',
      '.sc-path{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:11px;word-break:break-all;color:var(--dsw-alias-label-primary)}',

      // 状态 / 错误 / 崩溃
      '.sc-top-err{flex:none;padding:6px 12px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sc-recipe-strip{display:flex;align-items:center;gap:8px;flex:none;padding:5px 12px;border-bottom:1px solid var(--dsw-alias-border-l1);background:color-mix(in srgb,var(--dsw-alias-brand-primary) 7%,var(--dsw-alias-bg-layer-1));font-size:11.5px;color:var(--dsw-alias-brand-primary)}',
      '.sc-err{display:flex;flex-direction:column;gap:5px;padding:7px 9px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-error-primary) 45%,transparent);border-radius:10px;background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 8%,var(--dsw-alias-bg-layer-1))}',
      '.sc-err-title{font-size:12px;font-weight:600;color:var(--dsw-alias-state-error-primary)}',
      '.sc-err-msg{font-size:11.5px;white-space:pre-wrap;word-break:break-word;color:var(--dsw-alias-label-primary)}',
      '.sc-err-actions{display:flex;gap:5px}',
      '.sc-loading{display:flex;align-items:center;gap:7px;font-size:12px;color:var(--dsw-alias-label-secondary)}',
      '.sc-spinner{display:inline-block;width:12px;height:12px;flex:none;border:2px solid var(--dsw-alias-border-l2);border-top-color:var(--dsw-alias-brand-primary);border-radius:50%;animation:sc-spin .8s linear infinite}',
      '@keyframes sc-spin{to{transform:rotate(360deg)}}',
      '.sc-empty{display:flex;flex:1 1 auto;align-items:center;justify-content:center;padding:28px 20px;overflow-y:auto}',
      '.sc-empty-card{display:flex;flex-direction:column;gap:10px;width:100%;max-width:520px;padding:22px 20px;border:1px solid var(--dsw-alias-border-l1);border-radius:12px;background:var(--dsw-alias-bg-layer-1);text-align:center;align-items:center}',
      '.sc-empty-title{font-size:15px;font-weight:600}',
      '.sc-empty-hint{font-size:12px;color:var(--dsw-alias-label-secondary);max-width:420px}',
      '.sc-empty-actions{display:flex;gap:6px;flex-wrap:wrap;justify-content:center}',
      '.sc-crash{display:flex;flex-direction:column;align-items:flex-start;gap:7px;margin:14px;padding:16px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-error-primary) 45%,transparent);border-radius:12px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-primary)}',
      '.sc-crash-title{font-size:13.5px;font-weight:600;color:var(--dsw-alias-state-error-primary)}',
      '.sc-crash-hint{font-size:12px;color:var(--dsw-alias-label-secondary);max-width:520px}',
      '.sc-crash-detail{font-size:11px;color:var(--dsw-alias-label-secondary);word-break:break-word;max-height:120px;overflow:auto}',
      '.sc-crash-mini{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;border-radius:50%;border:1px solid var(--dsw-alias-state-error-primary);color:var(--dsw-alias-state-error-primary);font-size:11px;font-weight:700}',

      // v3：导入数据浮层
      '.sc-modal{position:fixed;inset:0;z-index:1300;display:flex;align-items:center;justify-content:center;padding:24px;background:color-mix(in srgb,var(--dsw-alias-label-primary) 34%,transparent)}',
      '.sc-modal-card{display:flex;flex-direction:column;gap:10px;width:min(780px,94vw);max-height:88vh;overflow-y:auto;padding:14px 16px 16px;border:1px solid var(--dsw-alias-border-l2);border-radius:12px;background:var(--dsw-alias-bg-layer-1);box-shadow:0 18px 48px color-mix(in srgb,var(--dsw-alias-label-primary) 24%,transparent)}',
      '.sc-modal-head{display:flex;align-items:flex-start;justify-content:space-between;gap:10px;flex:none}',
      '.sc-modal-title{font-size:14px;font-weight:600}',
      '.sc-modal-subtitle{font-size:11px;color:var(--dsw-alias-label-secondary);max-width:560px;margin-top:1px}',
      '.sc-modal-tabs{display:flex;gap:3px;flex:none;padding:2px;border:1px solid var(--dsw-alias-border-l1);border-radius:9px;background:var(--dsw-alias-bg-base);border-bottom:none}',
      '.sc-modal-tab{flex:1 1 0;height:26px;border:1px solid transparent;border-radius:7px;background:transparent;color:var(--dsw-alias-label-secondary);font:inherit;font-size:11.5px;cursor:pointer;white-space:nowrap}',
      '.sc-modal-tab:hover{background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary)}',
      '.sc-modal-tab-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary);font-weight:600}',
      '.sc-modal-section{display:flex;flex-direction:column;gap:5px;flex:none}',
      '.sc-modal-label{font-size:11.5px;color:var(--dsw-alias-label-secondary);font-weight:500}',
      '.sc-modal-hint{font-size:10.5px;color:var(--dsw-alias-label-tertiary);line-height:1.5}',
      '.sc-modal-input{width:100%;height:28px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);color:var(--dsw-alias-label-primary);font:inherit;font-size:12px;outline:none;min-width:0}',
      '.sc-modal-input:focus{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1))}',
      '.sc-form-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:8px;flex:none}',
      '.sc-form-cell{display:flex;flex-direction:column;gap:4px;min-width:0}',
      '.sc-file-row{display:flex;align-items:center;gap:8px;flex-wrap:wrap}',
      '.sc-file-chosen{font-size:11px;color:var(--dsw-alias-brand-primary);word-break:break-all}',
      '.sc-adv-toggle{align-self:flex-start}',
      '.sc-adv-body{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:8px;padding:8px 10px;border:1px dashed var(--dsw-alias-border-l1);border-radius:9px;background:var(--dsw-alias-bg-base);flex:none}',
      '.sc-cand-list{display:flex;flex-direction:column;gap:4px;max-height:220px;overflow-y:auto}',
      '.sc-cand{display:flex;align-items:center;justify-content:space-between;gap:8px;padding:5px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base)}',
      '.sc-cand-main{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.sc-cand-label{font-size:11.5px;font-weight:500;word-break:break-all}',
      '.sc-cand-detail{font-size:10.5px;color:var(--dsw-alias-label-tertiary);word-break:break-all}',
      '.sc-import-result{display:flex;flex-direction:column;gap:3px;padding:7px 9px;border:1px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 30%,var(--dsw-alias-border-l1));border-radius:9px;background:color-mix(in srgb,var(--dsw-alias-brand-primary) 5%,var(--dsw-alias-bg-base))}',
      '.sc-import-result-row{display:flex;align-items:baseline;gap:6px;flex-wrap:wrap;font-size:11px}',
      '.sc-import-result-k{color:var(--dsw-alias-label-tertiary)}',
      '.sc-import-result-v{color:var(--dsw-alias-label-primary);word-break:break-all}',
      '.sc-modal-actions{display:flex;align-items:center;justify-content:flex-end;gap:6px;flex-wrap:wrap;flex:none;padding-top:2px}',

      // v3：高级分析 tab
      '.sc-analysis{display:flex;flex:1 1 auto;gap:10px;min-height:0;min-width:0}',
      '.sc-analysis-list{display:flex;flex-direction:column;gap:4px;flex:0 0 236px;width:236px;min-height:0;padding:6px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1);overflow-y:auto}',
      '.sc-analysis-group{display:flex;flex-direction:column;gap:1px}',
      '.sc-analysis-group-title{font-size:10.5px;font-weight:600;letter-spacing:.03em;color:var(--dsw-alias-label-tertiary);padding:5px 7px 2px}',
      '.sc-analysis-item{display:flex;align-items:center;gap:6px;width:100%;padding:4px 7px;border:1px solid transparent;border-radius:7px;background:transparent;color:var(--dsw-alias-label-primary);font:inherit;font-size:11.5px;text-align:left;cursor:pointer;min-width:0}',
      '.sc-analysis-item:hover{background:var(--dsw-alias-bg-layer-2)}',
      '.sc-analysis-item-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 9%,var(--dsw-alias-bg-layer-1))}',
      '.sc-analysis-item-external{color:var(--dsw-alias-label-tertiary);cursor:help}',
      '.sc-analysis-item-title{flex:1 1 auto;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-analysis-tag{flex:none;font-size:9.5px;line-height:16px;padding:0 5px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;color:var(--dsw-alias-label-tertiary);white-space:nowrap}',
      '.sc-analysis-tag-native{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1));color:var(--dsw-alias-state-success-primary)}',
      '.sc-analysis-tag-external{border-style:dashed}',
      '.sc-analysis-tag-unavailable{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 40%,var(--dsw-alias-border-l1));color:var(--dsw-alias-state-warn-primary)}',
      '.sc-analysis-main{display:flex;flex:1 1 auto;flex-direction:column;gap:10px;min-width:0;min-height:0;overflow-y:auto}',
      '.sc-analysis-head{display:flex;flex-direction:column;gap:3px;flex:none}',
      '.sc-analysis-name{font-size:13.5px;font-weight:600}',
      '.sc-analysis-method{font-size:11px;color:var(--dsw-alias-label-secondary);word-break:break-word}',
      '.sc-analysis-needs{font-size:11px;color:var(--dsw-alias-state-warn-primary);word-break:break-word}',
      '.sc-analysis-note{font-size:11px;color:var(--dsw-alias-label-secondary);word-break:break-word}',
      '.sc-analysis-form{display:flex;flex-wrap:wrap;gap:8px;padding:9px 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.sc-analysis-field{display:flex;flex-direction:column;gap:3px;flex:1 1 170px;min-width:140px}',
      '.sc-analysis-field-label{font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-analysis-field-hint{font-size:10.5px;color:var(--dsw-alias-label-tertiary);line-height:1.45}',
      '.sc-analysis-checkbox{display:inline-flex;align-items:center;gap:5px;font-size:12px;color:var(--dsw-alias-label-primary);min-height:28px}',
      '.sc-analysis-actions{display:flex;align-items:center;gap:6px;flex-wrap:wrap;flex:1 1 100%}',
      '.sc-analysis-result{display:flex;flex-direction:column;gap:10px;padding:10px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.sc-analysis-block{display:flex;flex-direction:column;gap:5px;min-width:0}',
      '.sc-analysis-block-title{font-size:11.5px;font-weight:600;color:var(--dsw-alias-label-secondary)}',
      '.sc-analysis-bullets{margin:0;padding-left:18px;display:flex;flex-direction:column;gap:3px;font-size:12px;line-height:1.55}',
      '.sc-analysis-metrics{display:flex;flex-wrap:wrap;gap:5px}',
      '.sc-analysis-metric{display:inline-flex;align-items:baseline;gap:4px;padding:2px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:6px;background:var(--dsw-alias-bg-base);font-size:11px}',
      '.sc-analysis-metric-k{color:var(--dsw-alias-label-tertiary)}',
      '.sc-analysis-metric-v{font-weight:600}',
      '.sc-analysis-table-wrap{max-height:320px;overflow:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:8px}',
      '.sc-analysis-table{width:100%;border-collapse:collapse;font-size:11px}',
      '.sc-analysis-table th{position:sticky;top:0;z-index:1;padding:5px 8px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-secondary);font-weight:600;text-align:left;white-space:nowrap}',
      '.sc-analysis-table td{padding:3px 8px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 55%,transparent);white-space:nowrap;color:var(--dsw-alias-label-primary)}',
      '.sc-analysis-table tbody tr:last-child td{border-bottom:none}',
      '.sc-analysis-table-note{font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-analysis-figure{display:flex;flex-direction:column;gap:4px}',
      '.sc-analysis-figure img{display:block;max-width:100%;height:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-2)}',
      '.sc-analysis-figure-title{font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-analysis-empty{font-size:12px;color:var(--dsw-alias-label-secondary)}',
      '.sc-analysis-warn-list{display:flex;flex-direction:column;gap:2px;font-size:11px;color:var(--dsw-alias-state-warn-primary)}',

      // v4：样本 & 分组
      '.sc-merge-notice{display:flex;align-items:center;gap:8px;flex-wrap:wrap;flex:none;padding:6px 12px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-state-success-primary) 35%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 8%,var(--dsw-alias-bg-layer-1));font-size:11.5px;color:var(--dsw-alias-label-primary)}',
      '.sc-merge-notice-text{flex:1 1 auto;min-width:0}',
      '.sc-samples{display:flex;flex-direction:column;gap:6px}',
      '.sc-samples-head{display:flex;align-items:flex-start;justify-content:space-between;gap:8px;flex-wrap:wrap}',
      '.sc-samples-head-actions{display:flex;align-items:center;gap:5px;flex-wrap:wrap}',
      '.sc-samples-table-wrap{max-height:260px;overflow:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base)}',
      '.sc-samples-table{width:100%;border-collapse:collapse;font-size:11.5px}',
      '.sc-samples-table th{position:sticky;top:0;z-index:1;padding:5px 8px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-weight:600;text-align:left;white-space:nowrap}',
      '.sc-samples-table td{padding:3px 8px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 55%,transparent);vertical-align:middle}',
      '.sc-samples-table tbody tr:last-child td{border-bottom:none}',
      '.sc-samples-src{max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:var(--dsw-alias-label-secondary)}',
      '.sc-sample-name{width:150px;height:24px;font-size:11.5px}',
      '.sc-sample-group{width:120px;height:24px;font-size:11.5px}',
      '.sc-samples-status{font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-samples-tools{display:flex;flex-direction:column;gap:5px;padding-top:2px}',
      '.sc-samples-tool{display:flex;align-items:center;gap:6px;flex-wrap:wrap;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-samples-tool-label{flex:0 0 auto;min-width:120px;color:var(--dsw-alias-label-secondary)}',
      '.sc-samples-tool-input{width:250px;height:24px;font-size:11.5px}',
      '.sc-samples-key{width:140px}',
      '.sc-samples-preview{display:flex;flex-direction:column;gap:4px;padding:7px 9px;border:1px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 32%,var(--dsw-alias-border-l1));border-radius:8px;background:color-mix(in srgb,var(--dsw-alias-brand-primary) 5%,var(--dsw-alias-bg-base))}',
      '.sc-samples-preview-list{display:flex;flex-wrap:wrap;gap:4px;max-height:120px;overflow-y:auto}',
      '.sc-samples-preview-chip{display:inline-flex;align-items:center;gap:4px;padding:1px 7px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);font-size:10.5px}',
      '.sc-samples-preview-chip b{color:var(--dsw-alias-brand-primary);font-weight:600}',
      '.sc-samples-notice{display:flex;align-items:center;gap:8px;flex-wrap:wrap;padding:6px 9px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-success-primary) 35%,var(--dsw-alias-border-l1));border-radius:8px;background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 7%,var(--dsw-alias-bg-base));font-size:11.5px}',
      '.sc-samples-empty{font-size:11.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-path-lines{height:auto;min-height:64px;padding:6px 9px;line-height:1.5;font-size:11.5px;resize:vertical;font-family:inherit}',

      // v5：数据根目录扫描 → 样本表
      '.sc-import-step-title{font-size:12.5px;font-weight:600;color:var(--dsw-alias-label-primary);margin-bottom:2px}',
      '.sc-scan-row{display:flex;align-items:center;gap:6px;flex-wrap:wrap}',
      '.sc-scan-root{flex:1 1 260px;min-width:200px}',
      '.sc-scan-label{font-size:11.5px;color:var(--dsw-alias-label-secondary);white-space:nowrap}',
      '.sc-scan-toolbar{display:flex;align-items:center;gap:6px;flex-wrap:wrap;font-size:11px}',
      '.sc-scan-toolbar .sc-btn{height:24px;padding:0 8px;font-size:11px;border-radius:6px}',
      '.sc-scan-newgroup{width:150px;height:24px;font-size:11.5px}',
      '.sc-scan-newgroup-color{width:30px;height:24px;padding:1px;border:1px solid var(--dsw-alias-border-l1);border-radius:6px;background:var(--dsw-alias-bg-base)}',
      '.sc-scan-batch{display:flex;align-items:center;gap:5px}',
      '.sc-scan-batch select{height:24px;font-size:11.5px;max-width:150px}',
      '.sc-scan-groups{display:flex;flex-wrap:wrap;gap:4px}',
      '.sc-group-chip{display:inline-flex;align-items:center;padding:1px 8px;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);font-size:10.5px;color:var(--dsw-alias-label-primary);background:var(--dsw-alias-bg-layer-2)}',
      '.sc-group-chip-color{display:inline-block;width:7px;height:7px;border-radius:50%;margin-right:4px}',
      '.sc-scan-table-wrap{max-height:260px;overflow:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base)}',
      '.sc-scan-table{width:100%;border-collapse:collapse;font-size:11.5px}',
      '.sc-scan-table th{position:sticky;top:0;z-index:1;padding:5px 8px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-weight:600;text-align:left;white-space:nowrap}',
      '.sc-scan-table td{padding:3px 8px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 55%,transparent);vertical-align:middle}',
      '.sc-scan-table tbody tr:last-child td{border-bottom:none}',
      '.sc-scan-name{width:130px;height:24px;font-size:11.5px}',
      '.sc-scan-group-select{width:110px;height:24px;font-size:11.5px}',
      '.sc-scan-raw{max-width:170px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:var(--dsw-alias-label-secondary)}',
      '.sc-scan-detail{max-width:150px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:var(--dsw-alias-label-tertiary)}',
      '.sc-scan-geo{max-width:170px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:var(--dsw-alias-brand-primary)}',
      '.sc-scan-errors{display:flex;flex-direction:column;gap:2px;font-size:11px;color:var(--dsw-alias-state-warn-primary)}',
      '.sc-scan-unmatched{display:flex;flex-wrap:wrap;gap:4px;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-advanced-toggle{align-self:flex-start}',
      '.sc-advanced-section{display:flex;flex-direction:column;gap:8px;padding:8px 10px;border:1px dashed var(--dsw-alias-border-l2);border-radius:9px;background:var(--dsw-alias-bg-base)}',

      // v6：降维聚类控制卡 / 手动注释 / 分组维护 / Agent 输入
      '.sc-dimred{display:flex;flex-direction:column;gap:8px;padding:10px 12px;border:1px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 28%,var(--dsw-alias-border-l1));border-radius:10px;background:color-mix(in srgb,var(--dsw-alias-brand-primary) 4%,var(--dsw-alias-bg-layer-1));flex:none}',
      '.sc-dimred-head{display:flex;align-items:flex-start;justify-content:space-between;gap:8px;flex-wrap:wrap}',
      '.sc-dimred-title{font-size:12.5px;font-weight:600}',
      '.sc-dimred-sub{font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sc-dimred-grid{display:flex;flex-wrap:wrap;gap:8px 14px}',
      '.sc-dimred-field{display:flex;flex-direction:column;gap:3px;min-width:120px;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-dimred-field .sc-input{height:24px;font-size:11.5px;width:130px}',
      '.sc-dimred-slider-row{display:flex;align-items:center;gap:8px;flex:1 1 260px;min-width:220px;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sc-dimred-slider{flex:1 1 auto;min-width:120px;accent-color:var(--dsw-alias-brand-primary)}',
      '.sc-dimred-value{font-weight:600;color:var(--dsw-alias-brand-primary);min-width:34px}',
      '.sc-dimred-actions{display:flex;align-items:center;gap:6px;flex-wrap:wrap}',
      '.sc-dimred-warns{display:flex;flex-direction:column;gap:2px;font-size:10.5px;color:var(--dsw-alias-state-warn-primary)}',
      '.sc-annot{display:flex;flex-direction:column;gap:8px;padding:10px 12px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1);flex:none}',
      '.sc-annot-head{display:flex;align-items:flex-start;justify-content:space-between;gap:8px;flex-wrap:wrap}',
      '.sc-annot-title{font-size:12.5px;font-weight:600}',
      '.sc-annot-layout{display:flex;gap:10px;min-height:0;align-items:flex-start;flex-wrap:wrap}',
      '.sc-annot-clusters{display:flex;flex-direction:column;gap:3px;flex:0 0 210px;max-height:260px;overflow-y:auto;padding:5px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base)}',
      '.sc-annot-cluster{display:flex;flex-direction:column;gap:2px;padding:4px 6px;border:1px solid transparent;border-radius:6px;background:transparent;color:var(--dsw-alias-label-primary);font:inherit;font-size:11px;text-align:left;cursor:pointer}',
      '.sc-annot-cluster:hover{background:var(--dsw-alias-bg-layer-2)}',
      '.sc-annot-cluster-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 9%,var(--dsw-alias-bg-layer-1))}',
      '.sc-annot-cluster-top{color:var(--dsw-alias-label-tertiary);font-size:10px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sc-annot-main{display:flex;flex:1 1 320px;flex-direction:column;gap:6px;min-width:0}',
      '.sc-annot-genes{display:flex;flex-wrap:wrap;gap:4px}',
      '.sc-annot-gene-chip{display:inline-flex;align-items:center;padding:1px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary);font:inherit;font-size:10.5px;cursor:pointer}',
      '.sc-annot-gene-chip-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 55%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 12%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary);font-weight:600}',
      '.sc-annot-actions{display:flex;align-items:center;gap:5px;flex-wrap:wrap}',
      '.sc-annot-mapping{display:flex;flex-wrap:wrap;gap:4px;font-size:10.5px}',
      '.sc-annot-map-chip{padding:1px 8px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-success-primary) 35%,var(--dsw-alias-border-l1));border-radius:999px;background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 7%,var(--dsw-alias-bg-base))}',
      '.sc-annot-csv{display:flex;flex-direction:column;gap:4px}',
      '.sc-annot-csv textarea{width:100%;min-height:56px;padding:5px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);color:var(--dsw-alias-label-primary);font:inherit;font-size:11px;resize:vertical}',
      '.sc-groups-row{display:flex;flex-wrap:wrap;gap:4px;align-items:center}',
      '.sc-group-chip-del{margin-left:4px;border:none;background:transparent;color:var(--dsw-alias-state-error-primary);cursor:pointer;font-size:10.5px;padding:0 2px}',
      '.sc-group-chip-del:disabled{color:var(--dsw-alias-label-tertiary);cursor:default}',
      '.sc-sample-del{height:22px;padding:0 6px;font-size:10.5px;border-radius:6px}',
      '.sc-agent-ask{display:flex;flex-direction:column;gap:5px;padding:8px 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-base)}',
      '.sc-agent-ask textarea{width:100%;min-height:46px;padding:6px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-primary);font:inherit;font-size:11.5px;line-height:1.5;resize:vertical}',
      '.sc-agent-ask-actions{display:flex;align-items:center;justify-content:space-between;gap:6px}',
      '.sc-agent-sug{display:flex;flex-direction:column;gap:4px;padding:8px 10px;border:1px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 26%,var(--dsw-alias-border-l1));border-radius:10px;background:color-mix(in srgb,var(--dsw-alias-brand-primary) 4%,var(--dsw-alias-bg-base))}',
      '.sc-agent-sug-actions{display:flex;gap:5px;flex-wrap:wrap}',
      '.sc-agent-result{margin-top:2px;padding:4px 7px;border-radius:6px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-success-primary) 30%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 6%,transparent);font-size:10.5px;word-break:break-word}',

      '@media (max-width:1180px){.sc-left{flex-basis:240px;width:240px}.sc-right{flex-basis:262px;width:262px}}',
      '@media (max-width:980px){.sc-left{flex-basis:216px;width:216px}.sc-right{flex-basis:236px;width:236px}.sc-path-input{width:150px}}',
      '@media (prefers-reduced-motion:reduce){.sc-spinner{animation:none}}',
    ].join('')

    // ==========================================================================
    // 通用小组件
    // ==========================================================================

    function Spinner() {
      return h('span', { className: 'sc-spinner', 'aria-hidden': 'true' })
    }

    function LoadingLine(props) {
      return h('div', { className: 'sc-loading' }, h(Spinner), h('span', null, props.text || T('loading')))
    }

    function ErrorBanner(props) {
      return h('div', { className: 'sc-err', role: 'alert' },
        h('div', { className: 'sc-err-title' }, props.title || T('errorTitle')),
        props.message ? h('div', { className: 'sc-err-msg' }, String(props.message)) : null,
        props.onRetry
          ? h('div', { className: 'sc-err-actions' },
              h('button', { type: 'button', className: 'sc-btn sc-btn-sm', onClick: props.onRetry }, props.retryLabel || T('retry')))
          : null)
    }

    function Notice(props) {
      return h('div', { className: 'sc-notice' }, props.text || '')
    }

    /** 渲染崩溃兜底：每个主要区块都包一层，避免整个页面白屏。 */
    class PluginErrorBoundary extends React.Component {
      constructor(props) {
        super(props)
        this.state = { error: null }
        this.handleReset = this.handleReset.bind(this)
      }

      static getDerivedStateFromError(error) {
        return { error: error || new Error('unknown render error') }
      }

      componentDidCatch(error, info) {
        warnPlugin('UI crashed: ' + formatError(error), info)
      }

      handleReset() {
        this.setState({ error: null })
      }

      render() {
        if (this.state.error) {
          if (this.props.mini) return h('span', { className: 'sc-crash-mini', title: T('crashTitle') }, '!')
          return h('div', { className: 'sc-crash', role: 'alert' },
            h('div', { className: 'sc-crash-title' }, T('crashTitle')),
            h('div', { className: 'sc-crash-hint' }, T('crashHint')),
            h('div', { className: 'sc-crash-detail' }, formatError(this.state.error)),
            h('button', { type: 'button', className: 'sc-btn sc-btn-primary', onClick: this.handleReset }, T('retry')))
        }
        return this.props.children
      }
    }

    // ==========================================================================
    // 侧边栏图标（纯 SVG，不依赖外部资源）
    // ==========================================================================

    function SidebarGlyph(props) {
      var size = (props && props.size) || 18
      return h('svg', {
        width: size,
        height: size,
        viewBox: '0 0 24 24',
        fill: 'none',
        stroke: 'currentColor',
        strokeWidth: 1.6,
        strokeLinecap: 'round',
        strokeLinejoin: 'round',
        'aria-hidden': 'true',
        focusable: 'false',
      },
        h('circle', { cx: 9.2, cy: 9, r: 4 }),
        h('circle', { cx: 9.2, cy: 9, r: 1.2 }),
        h('circle', { cx: 16.8, cy: 15.8, r: 3 }),
        h('path', { d: 'M4 20.2h16' }),
        h('path', { d: 'M14 7.4l2.6 2.6' }))
    }

    function SidebarIcon(props) {
      return h(PluginErrorBoundary, { mini: true }, h(SidebarGlyph, { size: props && props.size }))
    }

    // ==========================================================================
    // 顶栏
    // ==========================================================================

    function EngineBadge(props) {
      var badge = engineBadge(props.engine)
      var label = T(badge.labelKey)
      var detail = badge.detail || ''
      var title = label + (detail ? ' · ' + detail : '')
      return h('span', { className: 'sc-badge sc-badge-' + badge.tone, title: title },
        h('span', { className: 'sc-badge-dot' }),
        label)
    }

    function TopBar(props) {
      var busy = props.busy || {}

      var datasetChips = (props.datasets || []).map(function (dataset) {
        return h('button', {
          key: dataset.id,
          type: 'button',
          className: 'sc-chip' + (dataset.id === props.active ? ' sc-chip-on' : ''),
          title: (dataset.name || dataset.id) + (dataset.sourceFormat ? ' · ' + dataset.sourceFormat : ''),
          disabled: Boolean(busy.select),
          onClick: function () { props.onSelect(dataset.id) },
        }, (dataset.name || dataset.id), dataset.nObs ? ' · ' + formatNumber(dataset.nObs) : '')
      })

      var recipeOptions = (props.recipes || []).map(function (recipe) {
        return h('option', { key: recipe.id, value: recipe.id }, recipe.title || recipe.id)
      })

      return h('header', { className: 'sc-top' },
        h('div', { className: 'sc-top-main' },
          h('div', { className: 'sc-title-block' },
            h('div', { className: 'sc-title' }, T('title')),
            h('div', { className: 'sc-subtitle' }, T('subtitle'))),
          h(EngineBadge, { engine: props.engine }),
          datasetChips.length
            ? h('div', { className: 'sc-chips', role: 'group', 'aria-label': T('title') }, datasetChips)
            : null),
        h('div', { className: 'sc-top-actions' },
          (props.recipes || []).length
            ? h('div', { className: 'sc-recipe' },
                h('select', {
                  className: 'sc-input sc-select sc-recipe-select',
                  value: props.recipeId || '',
                  disabled: props.recipeRunning,
                  'aria-label': T('recipes'),
                  title: T('recipes'),
                  onChange: function (event) { props.onRecipeChange(event.target.value) },
                },
                  h('option', { value: '' }, T('recipePick')),
                  recipeOptions),
                h('button', {
                  type: 'button',
                  className: 'sc-btn sc-btn-primary',
                  disabled: props.recipeRunning || !props.recipeId || !props.active,
                  onClick: props.onRunRecipe,
                }, props.recipeRunning ? T('recipeRunning') : T('runRecipe')))
            : null,
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-primary',
            disabled: Boolean(busy.import) || Boolean(busy.upload),
            title: T('importSubtitle'),
            onClick: props.onOpenImport,
          }, busy.import || busy.upload ? T('running') : T('importData')),
          h('button', {
            type: 'button',
            className: 'sc-btn',
            disabled: Boolean(busy.demo),
            onClick: props.onDemo,
          }, busy.demo ? T('running') : T('loadDemo')),
          h('button', {
            type: 'button',
            className: 'sc-btn',
            disabled: Boolean(busy.refresh),
            onClick: props.onRefresh,
          }, busy.refresh ? T('loading') : T('refresh'))))
    }


    // ==========================================================================
    // v3：导入数据面板（fixed 浮层；不走 react-dom portal，宿主只允许 require('react')）
    // ==========================================================================

    function ImportDialog(props) {
      // ---- 扫描流程（第 1 步）----
      var dirRef = React.useRef(null)
      var elapsedState = React.useState(0)
      var elapsed = elapsedState[0]
      var setElapsed = elapsedState[1]
      var startRef = React.useRef(0)
      var busy = Boolean(props.busy || props.inspectBusy || props.scanLoading || props.geoLoading)

      // ---- 高级（原有直接打开）流程本地状态 ----
      var advancedState = React.useState(false)
      var advanced = advancedState[0]
      var setAdvanced = advancedState[1]
      var sourceTypeState = React.useState('path')
      var sourceType = sourceTypeState[0]
      var setSourceType = sourceTypeState[1]
      var pathState = React.useState('')
      var path = pathState[0]
      var setPath = pathState[1]
      var mergeDirState = React.useState(false)
      var mergeDir = mergeDirState[0]
      var setMergeDir = mergeDirState[1]
      var urlState = React.useState('')
      var url = urlState[0]
      var setUrl = urlState[1]
      var fileState = React.useState(null)
      var file = fileState[0]
      var setFile = fileState[1]
      var formatState = React.useState('auto')
      var format = formatState[0]
      var setFormat = formatState[1]
      var transposeState = React.useState('auto')
      var transpose = transposeState[0]
      var setTranspose = transposeState[1]
      var objectState = React.useState('')
      var objectName = objectState[0]
      var setObjectName = objectState[1]
      var memberState = React.useState('')
      var member = memberState[0]
      var setMember = memberState[1]
      var fileRef = React.useRef(null)

      React.useEffect(function () {
        if (!busy) {
          setElapsed(0)
          return undefined
        }
        startRef.current = Date.now()
        var timer = setInterval(function () {
          setElapsed(Date.now() - startRef.current)
        }, 500)
        return function () { clearInterval(timer) }
      }, [busy])

      React.useEffect(function () {
        function onKeyDown(event) {
          if (event.key === 'Escape' && !busy && typeof props.onClose === 'function') props.onClose()
        }
        if (typeof document !== 'undefined' && document.addEventListener) document.addEventListener('keydown', onKeyDown)
        return function () {
          if (typeof document !== 'undefined' && document.removeEventListener) document.removeEventListener('keydown', onKeyDown)
        }
      }, [busy])

      // ---- 扫描流程动作 ----
      function pickDir() {
        if (busy) return
        // 优先让 Host 弹本机原生选择框（能拿到绝对路径）；老宿主才退回浏览器目录选择
        if (typeof props.onPickPath === 'function') {
          props.onPickPath()
          return
        }
        if (dirRef.current && dirRef.current.click) dirRef.current.click()
      }

      function handleDir(event) {
        var files = event && event.target ? event.target.files : null
        var list = files ? Array.prototype.slice.call(files) : []
        if (!list.length) return
        var first = list[0] || {}
        var full = typeof first.path === 'string' ? first.path : ''
        if (full) {
          var normalized = full.replace(/[\\/][^\\/]*$/, '')
          if (normalized) props.onScanRootChange(normalized)
        } else {
          props.onScanError(T('scanBrowseUnavailable'))
        }
        try {
          if (event && event.target) event.target.value = ''
        } catch (error) {
          /* ignore */
        }
      }

      function submitScan() {
        if (busy) return
        props.onScan()
      }

      var scanRows = props.scanRows || EMPTY_LIST
      var groups = props.scanGroups || EMPTY_LIST
      var groupColors = props.scanGroupColors || {}
      var scanResult = props.scanResult
      var geoResult = props.geoResult
      var geoQueried = Boolean(geoResult && geoResult.samples && geoResult.samples.length)
      var duplicates = findDuplicateSampleNames(scanRows)
      var selectedCount = 0
      for (var selectedIndex = 0; selectedIndex < scanRows.length; selectedIndex++) {
        if (scanRows[selectedIndex].selected === true) selectedCount += 1
      }
      var busyText = props.inspectBusy
        ? tf(T, 'importInspectBusy', formatDuration(elapsed) || '0 ms')
        : tf(T, 'importBusy', formatDuration(elapsed) || '0 ms')

      var scanBody = h(Fragment, null,
        h('div', { className: 'sc-import-step-title' }, T('scanStepTitle')),
        h('div', { className: 'sc-scan-row' },
          h('input', {
            className: 'sc-input sc-scan-root',
            type: 'text',
            value: props.scanRoot || '',
            placeholder: T('scanRootPlaceholder'),
            spellCheck: false,
            disabled: busy,
            'aria-label': T('scanRoot'),
            onChange: function (event) { props.onScanRootChange(event.target.value) },
            onKeyDown: function (event) {
              if (event.key === 'Enter') {
                event.preventDefault()
                submitScan()
              }
            },
          }),
          h('button', { type: 'button', className: 'sc-btn', disabled: busy, onClick: pickDir }, T('scanBrowse')),
          h('label', { className: 'sc-check' },
            h('input', {
              type: 'checkbox',
              checked: props.scanRecursive !== false,
              disabled: busy,
              onChange: function (event) { props.onScanRecursive(event.target.checked) },
            }),
            T('scanRecursive')),
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-primary',
            disabled: busy || !String(props.scanRoot || '').trim(),
            onClick: submitScan,
          }, props.scanLoading ? T('scanRunning') : T('scanRun'))),
        h('input', {
          ref: dirRef,
          type: 'file',
          multiple: true,
          webkitdirectory: '',
          directory: '',
          style: { display: 'none' },
          tabIndex: -1,
          'aria-hidden': 'true',
          onChange: handleDir,
        }),
        props.scanError ? h('div', { className: 'sc-inline-err' }, props.scanError) : null,
        scanResult && scanResult.errors && scanResult.errors.length
          ? h('div', { className: 'sc-scan-errors' }, scanResult.errors.slice(0, 8).map(function (item, index) {
              return h('div', { key: 'scan-err-' + index }, item)
            }))
          : null,
        h('div', { className: 'sc-scan-row' },
          h('span', { className: 'sc-scan-label' }, T('scanGeoTitle')),
          h('input', {
            className: 'sc-input sc-scan-newgroup',
            type: 'text',
            value: props.geoGse || '',
            placeholder: T('scanGeoPlaceholder'),
            spellCheck: false,
            disabled: busy,
            onChange: function (event) { props.onGeoGseChange(event.target.value) },
            onKeyDown: function (event) {
              if (event.key === 'Enter') {
                event.preventDefault()
                if (String(props.geoGse || '').trim()) props.onGeoQuery()
              }
            },
          }),
          h('button', {
            type: 'button',
            className: 'sc-btn',
            disabled: busy || !String(props.geoGse || '').trim(),
            onClick: props.onGeoQuery,
          }, props.geoLoading ? T('scanGeoRunning') : T('scanGeoRun'))),
        props.geoError ? h('div', { className: 'sc-inline-err' }, props.geoError) : null,
        geoResult
          ? h('div', { className: 'sc-samples-preview' },
              geoResult.gse
                ? h('div', { className: 'sc-scan-label' }, geoResult.gse + (geoResult.title ? ' · ' + geoResult.title : ''))
                : null,
              geoResult.samples && geoResult.samples.length
                ? h('div', { className: 'sc-hint-line' }, tf(T, 'scanGeoMatched', formatNumber(geoResult.matchedCount || 0)))
                : h('div', { className: 'sc-hint-line' }, T('scanGeoEmpty')),
              geoResult.unmatched && geoResult.unmatched.length
                ? h('div', { className: 'sc-scan-unmatched' },
                    T('scanGeoUnmatched') + '：' +
                    geoResult.unmatched.map(function (item) {
                      return item.gsm + (item.title ? '（' + item.title + '）' : '')
                    }).join('；'))
                : null)
          : null,
        scanRows.length
          ? h(Fragment, null,
              h('div', { className: 'sc-scan-toolbar' },
                h('button', { type: 'button', className: 'sc-btn', disabled: busy, onClick: props.onSelectAll }, T('scanSelectAll')),
                h('button', { type: 'button', className: 'sc-btn', disabled: busy, onClick: props.onClearAll }, T('scanClearAll')),
                h('span', { className: 'sc-muted' }, tf(T, 'scanSelectedCount', selectedCount + '/' + scanRows.length)),
                h('span', { className: 'sc-scan-label' }, T('scanNewGroup')),
                h('input', {
                  className: 'sc-input sc-scan-newgroup',
                  type: 'text',
                  value: props.scanNewGroup || '',
                  placeholder: T('scanNewGroupPlaceholder'),
                  disabled: busy,
                  onChange: function (event) { props.onScanNewGroupChange(event.target.value) },
                  onKeyDown: function (event) {
                    if (event.key === 'Enter') {
                      event.preventDefault()
                      props.onAddScanGroup()
                    }
                  },
                }),
                h('input', {
                  type: 'color',
                  className: 'sc-scan-newgroup-color',
                  value: props.scanNewGroupColor || '#4e79a7',
                  disabled: busy,
                  title: T('scanGroupColor'),
                  onChange: function (event) { props.onScanNewGroupColorChange(event.target.value) },
                }),
                h('button', {
                  type: 'button',
                  className: 'sc-btn',
                  disabled: busy || !String(props.scanNewGroup || '').trim(),
                  onClick: props.onAddScanGroup,
                }, T('scanAddGroup')),
                h('span', { className: 'sc-scan-label' }, T('scanBatchGroup')),
                h('select', {
                  className: 'sc-input sc-select',
                  value: props.scanBatchGroup || '',
                  disabled: busy,
                  onChange: function (event) { props.onBatchGroupChange(event.target.value) },
                },
                  h('option', { value: '' }, T('scanBatchPlaceholder')),
                  groups.map(function (group) {
                    return h('option', { key: 'batch-' + group, value: group }, group)
                  })),
                h('button', {
                  type: 'button',
                  className: 'sc-btn',
                  disabled: busy || !props.scanBatchGroup || !selectedCount,
                  onClick: props.onApplyBatchGroup,
                }, T('scanApplyBatch'))),
              duplicates.length
                ? h('div', { className: 'sc-inline-err' }, tf(T, 'scanDuplicateNames', duplicates.join('、')))
                : null,
              groups.length
                ? h('div', { className: 'sc-scan-groups' }, groups.map(function (group) {
                    return h('span', { key: 'chip-' + group, className: 'sc-group-chip' },
                      h('span', {
                        className: 'sc-group-chip-color',
                        style: { background: groupColors[group] || 'var(--dsw-alias-brand-primary)' },
                      }),
                      group)
                  }))
                : h('div', { className: 'sc-hint-line' }, T('scanGroupsEmpty')),
              h('div', { className: 'sc-scan-table-wrap' },
                h('table', { className: 'sc-scan-table' },
                  h('thead', null, h('tr', null,
                    h('th', null, T('scanColSelect')),
                    h('th', null, T('scanColName')),
                    h('th', null, T('scanColRaw')),
                    h('th', null, T('scanColFormat')),
                    h('th', null, T('scanColDetail')),
                    h('th', null, T('scanColGroup')),
                    geoQueried ? h('th', null, T('scanColGeo')) : null)),
                  h('tbody', null, scanRows.map(function (row) {
                    var optionNodes = [h('option', { key: 'none', value: '' }, T('scanGroupPlaceholder'))]
                    for (var groupIndex = 0; groupIndex < groups.length; groupIndex++) {
                      optionNodes.push(h('option', { key: 'g-' + groups[groupIndex], value: groups[groupIndex] }, groups[groupIndex]))
                    }
                    if (row.group && groups.indexOf(row.group) < 0) {
                      optionNodes.push(h('option', { key: 'g-current', value: row.group }, row.group))
                    }
                    return h('tr', { key: 'scan-' + row.key },
                      h('td', null, h('input', {
                        type: 'checkbox',
                        checked: row.selected === true,
                        disabled: busy,
                        'aria-label': T('scanColSelect'),
                        onChange: function (event) { props.onScanRowSelect(row.key, event.target.checked) },
                      })),
                      h('td', null, h('input', {
                        className: 'sc-input sc-scan-name',
                        type: 'text',
                        value: row.name || '',
                        disabled: busy,
                        'aria-label': T('scanColName'),
                        onChange: function (event) { props.onScanRowName(row.key, event.target.value) },
                      })),
                      h('td', { className: 'sc-scan-raw', title: row.rawName }, row.rawName),
                      h('td', null, row.formatLabel || row.format || '—'),
                      h('td', { className: 'sc-scan-detail', title: row.detail || row.path }, row.detail || row.path || '—'),
                      h('td', null, h('select', {
                        className: 'sc-input sc-select sc-scan-group-select',
                        value: row.group || '',
                        disabled: busy,
                        'aria-label': T('scanColGroup'),
                        onChange: function (event) { props.onScanRowGroup(row.key, event.target.value) },
                      }, optionNodes)),
                      geoQueried
                        ? h('td', { className: 'sc-scan-geo', title: row.geoTitle || '' }, row.geoTitle || '—')
                        : null)
                  })))))
          : h('div', { className: 'sc-notice' }, T('scanEmpty')),
        props.scanConfirmError ? h('div', { className: 'sc-inline-err' }, props.scanConfirmError) : null,
        h('div', { className: 'sc-modal-actions' },
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-primary',
            disabled: busy || !selectedCount,
            onClick: props.onConfirmSamples,
          }, props.busy ? T('scanConfirmBusy') : T('scanConfirm'))))

      // ---- 高级：原有直接打开文件/URL（单文件、多路径、目录合并、压缩包、RData）----
      function buildOptions() {
        var lines = sourceType === 'path' ? parsePathLines(path) : []
        return {
          sourceType: sourceType,
          path: sourceType === 'path' && lines.length <= 1 ? (lines[0] || '') : '',
          paths: sourceType === 'path' ? lines : EMPTY_LIST,
          mergeDir: sourceType === 'path' ? mergeDir : false,
          url: sourceType === 'url' ? url : '',
          file: sourceType === 'file' ? file : null,
          format: format,
          transpose: transpose,
          object: objectName,
          member: member,
        }
      }

      function canSubmit() {
        if (busy) return false
        if (sourceType === 'path') return parsePathLines(path).length > 0
        if (sourceType === 'url') return Boolean(String(url || '').trim())
        return Boolean(file)
      }

      function submitAdvanced(event) {
        if (event && event.preventDefault) event.preventDefault()
        if (!canSubmit()) {
          props.onLocalError(T('importNeedSource'))
          return
        }
        props.onImport(buildOptions())
      }

      function inspectAdvanced() {
        if (busy) return
        if (sourceType === 'file') {
          props.onLocalError(T('importInspectFileHint'))
          return
        }
        var options = buildOptions()
        var inspectPath = options.path || (options.paths && options.paths.length ? options.paths[0] : '')
        if (!inspectPath && !options.url.trim()) {
          props.onLocalError(T('importNeedSource'))
          return
        }
        props.onInspect({
          sourceType: options.sourceType,
          path: inspectPath,
          paths: options.paths,
          url: options.url,
          format: options.format,
          member: options.member,
          object: options.object,
          transpose: options.transpose,
          mergeDir: options.mergeDir,
        })
      }

      function pickFile() {
        if (busy) return
        if (fileRef.current && fileRef.current.click) fileRef.current.click()
      }

      function handleFile(event) {
        var files = event && event.target ? event.target.files : null
        var chosen = files && files.length ? files[0] : null
        if (chosen) {
          setFile(chosen)
          setSourceType('file')
        }
        try {
          if (event && event.target) event.target.value = ''
        } catch (error) {
          /* ignore */
        }
      }

      var formats = props.formats && props.formats.length ? props.formats : normalizeFormats(null)
      var formatOptions = formats.map(function (item) {
        return h('option', {
          key: item.id,
          value: item.id,
          disabled: item.available === false,
          title: item.note || '',
        }, item.title + (item.available === false ? ' · ' + T('analysisUnavailable') : ''))
      })
      var candidates = props.candidates || EMPTY_LIST
      var inspectResult = props.inspectResult
      var sourceTabs = [
        { id: 'path', labelKey: 'importSourcePath' },
        { id: 'url', labelKey: 'importSourceUrl' },
        { id: 'file', labelKey: 'importSourceFile' },
      ].map(function (tab) {
        return h('button', {
          key: tab.id,
          type: 'button',
          disabled: busy,
          className: 'sc-modal-tab' + (sourceType === tab.id ? ' sc-modal-tab-on' : ''),
          onClick: function () { setSourceType(tab.id) },
        }, T(tab.labelKey))
      })

      var sourceBody = null
      if (sourceType === 'path') {
        sourceBody = h('div', { className: 'sc-modal-section' },
          h('label', { className: 'sc-modal-label', htmlFor: 'sc-import-path' }, T('importPathLines')),
          h('textarea', {
            id: 'sc-import-path',
            className: 'sc-modal-input sc-path-lines',
            rows: 3,
            value: path,
            placeholder: T('importPathLinesPlaceholder'),
            spellCheck: false,
            disabled: busy,
            onChange: function (event) { setPath(event.target.value) },
          }),
          h('label', { className: 'sc-check' },
            h('input', {
              type: 'checkbox',
              checked: mergeDir,
              disabled: busy,
              onChange: function (event) { setMergeDir(event.target.checked) },
            }),
            T('importMergeDir')),
          h('div', { className: 'sc-modal-hint' }, T('importMultiHint')),
          h('div', { className: 'sc-modal-hint' }, T('importPathHint')))
      } else if (sourceType === 'url') {
        sourceBody = h('div', { className: 'sc-modal-section' },
          h('label', { className: 'sc-modal-label', htmlFor: 'sc-import-url' }, T('importSourceUrl')),
          h('input', {
            id: 'sc-import-url',
            className: 'sc-modal-input',
            type: 'text',
            value: url,
            placeholder: T('importUrlPlaceholder'),
            spellCheck: false,
            disabled: busy,
            onChange: function (event) { setUrl(event.target.value) },
          }),
          h('div', { className: 'sc-modal-hint' }, T('importUrlHint')))
      } else {
        sourceBody = h('div', { className: 'sc-modal-section' },
          h('div', { className: 'sc-file-row' },
            h('button', { type: 'button', className: 'sc-btn', disabled: busy, onClick: pickFile }, T('importChooseFile')),
            file ? h('span', { className: 'sc-file-chosen' }, tf(T, 'importFileChosen', file.name || '')) : null),
          h('div', { className: 'sc-modal-hint' }, T('importFileUploadHint')))
      }

      var optionsState = React.useState(false)
      var optionsOpen = optionsState[0]
      var setOptionsOpen = optionsState[1]
      var advancedOptions = null
      if (optionsOpen) {
        advancedOptions = h('div', { className: 'sc-adv-body' },
          h('div', { className: 'sc-form-cell' },
            h('label', { className: 'sc-modal-label' }, T('importTranspose')),
            h('select', {
              className: 'sc-input sc-select sc-modal-input',
              value: transpose,
              disabled: busy,
              onChange: function (event) { setTranspose(event.target.value) },
            },
              h('option', { value: 'auto' }, T('importTransposeAuto')),
              h('option', { value: 'false' }, T('importTransposeNo')),
              h('option', { value: 'true' }, T('importTransposeYes')))),
          h('div', { className: 'sc-form-cell' },
            h('label', { className: 'sc-modal-label' }, T('importObject')),
            h('input', {
              className: 'sc-modal-input',
              type: 'text',
              value: objectName,
              placeholder: T('importObjectPlaceholder'),
              disabled: busy,
              spellCheck: false,
              onChange: function (event) { setObjectName(event.target.value) },
            })),
          h('div', { className: 'sc-form-cell' },
            h('label', { className: 'sc-modal-label' }, T('importMember')),
            h('input', {
              className: 'sc-modal-input',
              type: 'text',
              value: member,
              placeholder: T('importMemberPlaceholder'),
              disabled: busy,
              spellCheck: false,
              onChange: function (event) { setMember(event.target.value) },
            })))
      }

      var candidateNodes = null
      if (candidates.length) {
        candidateNodes = h('div', { className: 'sc-modal-section' },
          h('div', { className: 'sc-modal-label' }, T('importCandidates')),
          props.candidatesDisabled ? h('div', { className: 'sc-modal-hint' }, T('importCandidateUploadHint')) : null,
          h('div', { className: 'sc-cand-list' }, candidates.map(function (candidate, index) {
            return h('div', { key: candidate.key || ('cand-' + index), className: 'sc-cand' },
              h('div', { className: 'sc-cand-main' },
                h('div', { className: 'sc-cand-label' }, candidate.label || '—'),
                candidate.detail ? h('div', { className: 'sc-cand-detail' }, candidate.detail) : null),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm sc-btn-primary',
                disabled: busy || Boolean(props.candidatesDisabled),
                onClick: function () { props.onPickCandidate(candidate) },
              }, T('importSelectCandidate')))
          })))
      }

      var inspectNodes = null
      if (inspectResult) {
        inspectNodes = h('div', { className: 'sc-modal-section' },
          h('div', { className: 'sc-modal-label' }, T('importInspectResult')),
          h('div', { className: 'sc-import-result' },
            h('div', { className: 'sc-import-result-row' },
              h('span', { className: 'sc-import-result-k' }, T('importInspectFormat')),
              h('span', { className: 'sc-import-result-v' }, inspectResult.format || '—')),
            h('div', { className: 'sc-import-result-row' },
              h('span', { className: 'sc-import-result-k' }, T('importInspectSize')),
              h('span', { className: 'sc-import-result-v' }, formatBytes(inspectResult.size))),
            h('div', { className: 'sc-import-result-row' },
              h('span', { className: 'sc-import-result-k' }, T('importInspectDims')),
              h('span', { className: 'sc-import-result-v' },
                (inspectResult.nObs ? formatNumber(inspectResult.nObs) : '?') + ' × ' + (inspectResult.nVars ? formatNumber(inspectResult.nVars) : '?'))),
            inspectResult.notes.length
              ? h('div', { className: 'sc-modal-hint' }, inspectResult.notes.join('；'))
              : null,
            inspectResult.candidates.length
              ? h('div', { className: 'sc-cand-list' }, inspectResult.candidates.map(function (candidate, index) {
                  return h('div', { key: 'inspect-cand-' + index, className: 'sc-cand' },
                    h('div', { className: 'sc-cand-main' },
                      h('div', { className: 'sc-cand-label' }, candidate.label || '—'),
                      candidate.detail ? h('div', { className: 'sc-cand-detail' }, candidate.detail) : null))
                }))
              : null))
      }

      var advancedBody = null
      if (advanced) {
        advancedBody = h('div', { className: 'sc-advanced-section' },
          h('div', { className: 'sc-modal-hint' }, T('scanAdvancedHint')),
          h('div', { className: 'sc-modal-tabs' }, sourceTabs),
          sourceBody,
          h('input', {
            ref: fileRef,
            type: 'file',
            accept: IMPORT_ACCEPT,
            style: { display: 'none' },
            tabIndex: -1,
            'aria-hidden': 'true',
            onChange: handleFile,
          }),
          h('div', { className: 'sc-form-grid' },
            h('div', { className: 'sc-form-cell' },
              h('label', { className: 'sc-modal-label' }, T('importFormat')),
              h('select', {
                className: 'sc-input sc-select sc-modal-input',
                value: format,
                disabled: busy,
                onChange: function (event) { setFormat(event.target.value) },
              }, formatOptions),
              props.formatsError
                ? h('div', { className: 'sc-modal-hint' },
                    T('importFormatsFallback'),
                    ' ',
                    h('button', {
                      type: 'button',
                      className: 'sc-btn sc-btn-sm',
                      disabled: busy,
                      onClick: props.onReloadFormats,
                    }, T('importRetryFormats')))
                : null)),
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-sm sc-adv-toggle',
            'aria-expanded': optionsOpen ? 'true' : 'false',
            onClick: function () { setOptionsOpen(!optionsOpen) },
          }, (optionsOpen ? '▾ ' : '▸ ') + T('importAdvanced')),
          advancedOptions,
          props.error ? h('div', { className: 'sc-inline-err' }, props.error) : null,
          candidateNodes,
          props.inspectError ? h('div', { className: 'sc-inline-err' }, props.inspectError) : null,
          inspectNodes,
          h('div', { className: 'sc-modal-actions' },
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-primary',
              disabled: !canSubmit(),
              onClick: submitAdvanced,
            }, props.busy ? busyText : T('importSubmit')),
            h('button', {
              type: 'button',
              className: 'sc-btn',
              disabled: busy || sourceType === 'file',
              onClick: inspectAdvanced,
            }, props.inspectBusy ? busyText : T('importInspect'))))
      }

      return h('div', {
        className: 'sc-modal',
        role: 'dialog',
        'aria-modal': 'true',
        'aria-label': T('importTitle'),
      },
        h('div', { className: 'sc-modal-card' },
          h('div', { className: 'sc-modal-head' },
            h('div', null,
              h('div', { className: 'sc-modal-title' }, T('importTitle')),
              h('div', { className: 'sc-modal-subtitle' }, T('importSubtitle'))),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: busy,
              onClick: props.onClose,
            }, T('importClose'))),
          scanBody,
          h('div', { className: 'sc-divider' }),
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-sm sc-advanced-toggle',
            'aria-expanded': advanced ? 'true' : 'false',
            onClick: function () { setAdvanced(!advanced) },
          }, (advanced ? '▾ ' : '▸ ') + T('scanAdvancedToggle')),
          advancedBody))
    }


    // ==========================================================================
    // v3：高级分析 tab
    // ==========================================================================

    function renderAnalysisParamControl(descriptor, value, onChange, obsColumns, disabled) {
      if (descriptor.type === 'number') {
        return h('input', {
          className: 'sc-input sc-modal-input',
          type: 'number',
          min: descriptor.min === null ? undefined : descriptor.min,
          max: descriptor.max === null ? undefined : descriptor.max,
          step: descriptor.step === null ? 'any' : descriptor.step,
          value: value === null || value === undefined ? '' : String(value),
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value) },
        })
      }
      if (descriptor.type === 'boolean') {
        // 外层 .sc-analysis-field 已有 label 文本，这里只给勾选框（点 label 也能切换）。
        return h('input', {
          className: 'sc-analysis-checkbox',
          type: 'checkbox',
          checked: value === true,
          disabled: disabled,
          onChange: function (event) { onChange(event.target.checked) },
        })
      }
      if (descriptor.type === 'select') {
        var optionValues = descriptor.options.slice()
        var current = value === null || value === undefined ? '' : String(value)
        var hasCurrent = false
        for (var index = 0; index < optionValues.length; index++) {
          if (optionValues[index].value === current) hasCurrent = true
        }
        if (current && !hasCurrent) optionValues.unshift({ value: current, label: current })
        return h('select', {
          className: 'sc-input sc-select sc-modal-input',
          value: current,
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value) },
        }, optionValues.map(function (option) {
          return h('option', { key: option.value, value: option.value }, option.label)
        }))
      }
      if (descriptor.type === 'column') {
        var columns = obsColumns || []
        var columnOptions = [{ value: '', label: T('analysisNoSelection') }]
        for (var columnIndex = 0; columnIndex < columns.length; columnIndex++) {
          var column = columns[columnIndex]
          columnOptions.push({
            value: column.key,
            label: column.key + (column.kind === 'numeric' ? ' · numeric' : ' · categorical'),
          })
        }
        var columnValue = value === null || value === undefined ? '' : String(value)
        var known = false
        for (var optionIndex = 0; optionIndex < columnOptions.length; optionIndex++) {
          if (columnOptions[optionIndex].value === columnValue) known = true
        }
        if (columnValue && !known) columnOptions.push({ value: columnValue, label: columnValue })
        return h('select', {
          className: 'sc-input sc-select sc-modal-input',
          value: columnValue,
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value) },
        }, columnOptions.map(function (option) {
          return h('option', { key: option.value, value: option.value }, option.label)
        }))
      }
      return h('input', {
        className: 'sc-input sc-modal-input',
        type: 'text',
        value: value === null || value === undefined ? '' : String(value),
        placeholder: descriptor.type === 'file' ? '/absolute/path/file.txt' : '',
        disabled: disabled,
        spellCheck: false,
        onChange: function (event) { onChange(event.target.value) },
      })
    }

    function renderAnalysisTableNode(table, index) {
      var view = formatTableRows(table, 200)
      return h('div', { key: 'table-' + (view.id || index), className: 'sc-analysis-block' },
        view.title ? h('div', { className: 'sc-analysis-figure-title' }, view.title) : null,
        h('div', { className: 'sc-analysis-table-wrap' },
          h('table', { className: 'sc-analysis-table' },
            h('thead', null, h('tr', null, view.columns.map(function (column) {
              return h('th', { key: column.key, className: column.numeric ? 'sc-td-num' : '' }, column.label)
            }))),
            h('tbody', null, view.rows.map(function (row) {
              return h('tr', { key: row.key }, row.cells.map(function (cell) {
                return h('td', { key: cell.key, className: cell.numeric ? 'sc-td-num' : '' }, cell.text)
              }))
            })))),
        h('div', { className: 'sc-analysis-table-note' },
          tf(T, 'analysisTableRows', view.total) +
          (view.truncated ? ' · ' + tf(T, 'analysisTableShown', view.shown) : '') +
          (view.note ? ' · ' + view.note : '')))
    }

    function renderAnalysisResultNode(result) {
      var blocks = []
      if (result.obsMerged.length) {
        blocks.push(h('div', { key: 'obs', className: 'sc-analysis-note' },
          tf(T, 'analysisNewColumns', result.obsMerged.join(', '))))
      }
      if (result.summary.bullets.length) {
        blocks.push(h('div', { key: 'bullets', className: 'sc-analysis-block' },
          h('div', { className: 'sc-analysis-block-title' }, T('analysisResult')),
          h('ul', { className: 'sc-analysis-bullets' }, result.summary.bullets.map(function (bullet, index) {
            return h('li', { key: 'b-' + index }, bullet)
          }))))
      }
      var metricKeys = Object.keys(result.summary.metrics)
      if (metricKeys.length) {
        blocks.push(h('div', { key: 'metrics', className: 'sc-analysis-block' },
          h('div', { className: 'sc-analysis-block-title' }, T('analysisMetrics')),
          h('div', { className: 'sc-analysis-metrics' }, metricKeys.map(function (key) {
            return h('span', { key: key, className: 'sc-analysis-metric' },
              h('span', { className: 'sc-analysis-metric-k' }, key),
              h('span', { className: 'sc-analysis-metric-v' }, formatTableCell(result.summary.metrics[key], typeof result.summary.metrics[key] === 'number')))
          }))))
      }
      if (result.tables.length) {
        blocks.push(h('div', { key: 'tables', className: 'sc-analysis-block' },
          h('div', { className: 'sc-analysis-block-title' }, T('analysisTables') + ' · ' + result.tables.length),
          result.tables.map(renderAnalysisTableNode)))
      }
      if (result.figures.length) {
        blocks.push(h('div', { key: 'figures', className: 'sc-analysis-block' },
          h('div', { className: 'sc-analysis-block-title' }, T('analysisFigures') + ' · ' + result.figures.length),
          result.figures.map(function (figure, index) {
            return h('div', { key: 'figure-' + (figure.id || index), className: 'sc-analysis-figure' },
              figure.title ? h('div', { className: 'sc-analysis-figure-title' }, figure.title) : null,
              figure.url
                ? h('img', { src: figure.url, alt: figure.title || T('analysisFigures'), loading: 'lazy' })
                : h(Fragment, null,
                    h('div', { className: 'sc-analysis-empty' }, T('analysisFigureNoPreview')),
                    h('div', { className: 'sc-path' }, figure.file || '—')))
          })))
      }
      if (result.warnings.length) {
        blocks.push(h('div', { key: 'warnings', className: 'sc-analysis-block' },
          h('div', { className: 'sc-analysis-block-title' }, T('analysisWarnings')),
          h('div', { className: 'sc-analysis-warn-list' }, result.warnings.map(function (warning, index) {
            return h('div', { key: 'w-' + index }, warning)
          }))))
      }
      if (!blocks.length) return h('div', { className: 'sc-analysis-empty' }, T('analysisEmptyResult'))
      return h('div', { className: 'sc-analysis-result' }, blocks)
    }

    function AnalysisTab(props) {
      if (!props.activeId) return Notice({ text: T('noData') })
      var groups = props.groups || EMPTY_LIST
      var selected = props.selected || null

      var listNode = null
      if (props.loading && !groups.length) {
        listNode = LoadingLine({ text: T('analysisLoading') })
      } else if (props.listError && !groups.length) {
        listNode = ErrorBanner({ message: props.listError, onRetry: props.onReload, retryLabel: T('retry') })
      } else if (!groups.length) {
        listNode = Notice({ text: T('analysisEmpty') })
      } else {
        listNode = h('div', { className: 'sc-analysis-list' }, groups.map(function (entry) {
          return h('div', { key: entry.group, className: 'sc-analysis-group' },
            h('div', { className: 'sc-analysis-group-title' }, entry.labelKey ? T(entry.labelKey) : entry.group),
            entry.items.map(function (item) {
              var tag = item.external
                ? { className: 'sc-analysis-tag sc-analysis-tag-external', text: T('analysisExternalShort') }
                : (item.available
                    ? { className: 'sc-analysis-tag sc-analysis-tag-native', text: T('analysisNativeShort') }
                    : { className: 'sc-analysis-tag sc-analysis-tag-unavailable', text: T('analysisUnavailable') })
              return h('button', {
                key: item.id,
                type: 'button',
                className: 'sc-analysis-item' +
                  (selected && selected.id === item.id ? ' sc-analysis-item-on' : '') +
                  (item.external ? ' sc-analysis-item-external' : ''),
                title: item.external ? (item.reason || item.title) : item.title,
                onClick: function () { props.onSelect(item) },
              },
                h('span', { className: 'sc-analysis-item-title' }, item.title),
                h('span', { className: tag.className }, tag.text))
            }))
        }))
      }

      var mainNode = null
      if (!selected) {
        mainNode = h('div', { className: 'sc-guide' },
          h('div', { className: 'sc-guide-hint' }, T('analysisSelectHint')))
      } else {
        var hints = needsHint(selected, { steps: props.doneSteps || EMPTY_LIST, obsColumns: props.obsColumns || EMPTY_LIST })
        var headNodes = [
          h('div', { key: 'name', className: 'sc-analysis-name' }, selected.title),
          selected.method
            ? h('div', { key: 'method', className: 'sc-analysis-method' }, T('analysisMethod') + ' · ' + selected.method)
            : null,
          selected.external
            ? h('div', { key: 'external', className: 'sc-analysis-needs' },
                T('analysisExternalReason') + '：' + (selected.reason || '—'))
            : null,
          selected.external && selected.alternative
            ? h('div', { key: 'alternative', className: 'sc-analysis-note' },
                T('analysisExternalAlternative') + '：' + selected.alternative)
            : null,
          !selected.external && !selected.available && selected.reason
            ? h('div', { key: 'unavailable', className: 'sc-analysis-needs' },
                T('analysisUnavailable') + '：' + selected.reason)
            : null,
          hints ? h('div', { key: 'needs', className: 'sc-analysis-needs' }, T('analysisNeeds') + '：' + hints) : null,
        ]
        var bodyNodes = [h('div', { key: 'head', className: 'sc-analysis-head' }, headNodes)]
        if (!selected.external) {
          var fields = (selected.params || []).map(function (descriptor) {
            return h('label', { key: descriptor.key, className: 'sc-analysis-field' },
              h('span', { className: 'sc-analysis-field-label' }, descriptor.label),
              renderAnalysisParamControl(
                descriptor,
                props.values ? props.values[descriptor.key] : undefined,
                function (value) { props.onValueChange(descriptor.key, value) },
                props.obsColumns || EMPTY_LIST,
                Boolean(props.running)),
              descriptor.hint ? h('span', { className: 'sc-analysis-field-hint' }, descriptor.hint) : null)
          })
          bodyNodes.push(h('div', { key: 'form', className: 'sc-analysis-form' },
            fields.length ? fields : h('div', { className: 'sc-analysis-empty' }, T('analysisNoParams')),
            h('div', { className: 'sc-analysis-actions' },
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-primary',
                disabled: Boolean(props.running) || (selected.available === false),
                onClick: props.onRun,
              }, props.running ? T('analysisRunning') : T('analysisRun')))))
        }
        if (props.error) {
          bodyNodes.push(h('div', { key: 'run-error' },
            ErrorBanner({ message: props.error, onRetry: props.onRun, retryLabel: T('retry') })))
        }
        if (props.result) bodyNodes.push(h('div', { key: 'result' }, renderAnalysisResultNode(props.result)))
        mainNode = h(Fragment, null, bodyNodes)
      }

      return h('div', { className: 'sc-analysis' },
        listNode,
        h('div', { className: 'sc-analysis-main' }, mainNode))
    }


    // ==========================================================================
    // 左栏：流程卡片
    // ==========================================================================

    function renderParamControl(step, param, value, onChange, disabled) {
      var inputId = 'sc-param-' + step.id + '-' + param.key
      if (param.kind === 'select') {
        var options = (param.options || []).map(function (option) {
          var optionValue = typeof option === 'string' ? option : option.value
          var optionLabel = typeof option === 'string'
            ? option
            : (option.labelKey ? T(option.labelKey) : String(option.value))
          return h('option', { key: optionValue, value: optionValue }, optionLabel)
        })
        return h('select', {
          id: inputId,
          className: 'sc-input sc-select',
          value: value === null || value === undefined ? '' : String(value),
          disabled: disabled,
          onChange: function (event) { onChange(param.key, event.target.value) },
        }, options)
      }
      if (param.kind === 'text') {
        return h('input', {
          id: inputId,
          className: 'sc-input',
          type: 'text',
          value: value === null || value === undefined ? '' : String(value),
          disabled: disabled,
          onChange: function (event) { onChange(param.key, event.target.value) },
        })
      }
      return h('input', {
        id: inputId,
        className: 'sc-input',
        type: 'number',
        min: param.min,
        max: param.max,
        step: param.step || 1,
        value: value === null || value === undefined ? '' : String(value),
        disabled: disabled,
        'aria-label': T(param.labelKey),
        onChange: function (event) { onChange(param.key, event.target.value) },
      })
    }

    function StepCard(props) {
      var expandedState = React.useState(false)
      var expanded = expandedState[0]
      var setExpanded = expandedState[1]
      var valuesState = React.useState(function () { return defaultParams(props.step.id) })
      var values = valuesState[0]
      var setValues = valuesState[1]
      var localErrorState = React.useState('')
      var localError = localErrorState[0]
      var setLocalError = localErrorState[1]

      var step = props.step
      var status = props.status || 'pending'
      var disabled = Boolean(props.disabled)
      var busy = Boolean(props.busy)

      function setValue(key, value) {
        setLocalError('')
        setValues(function (previous) {
          var next = Object.assign({}, previous)
          next[key] = value
          return next
        })
      }

      function toggleExpanded() {
        if (expanded) {
          setExpanded(false)
          return
        }
        setExpanded(true)
      }

      function handleKeyDown(event) {
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault()
          toggleExpanded()
        }
      }

      function handleRun() {
        if (busy || disabled) return
        if (step.id === 'annotate') {
          var mapping = props.annotationMapping || {}
          var keys = Object.keys(mapping)
          if (!keys.length) {
            setLocalError(T('annotationNeedMapping'))
            return
          }
          setLocalError('')
          props.onRun({ mapping: mapping })
          return
        }
        var built = collectStepParams(step.id, values)
        if (built.error) {
          setLocalError(built.error)
          return
        }
        setLocalError('')
        props.onRun(built.params)
      }

      var statusKey = stepStatusKey(status)
      var body = null
      if (expanded) {
        var rows = step.params.map(function (param) {
          return h('label', { key: param.key, className: 'sc-param' },
            h('span', { className: 'sc-param-label' }, T(param.labelKey)),
            renderParamControl(step, param, values[param.key], setValue, busy || disabled))
        })
        if (step.id === 'annotate') {
          rows.push(h('div', { key: 'annotate-hint', className: 'sc-param-hint' },
            tf(T, 'annotationCountHint', props.annotationCount || 0)))
        }
        if (step.requires.length) {
          rows.push(h('div', { key: 'requires', className: 'sc-param-hint' },
            T('needPrereq') + '：' + step.requires.map(stepTitleOf).join(' → ')))
        }
        body = h('div', { className: 'sc-step-params' }, rows)
      }

      return h('section', { className: 'sc-step sc-step-' + status },
        h('div', {
          className: 'sc-step-head',
          role: 'button',
          tabIndex: 0,
          'aria-expanded': expanded ? 'true' : 'false',
          onClick: toggleExpanded,
          onKeyDown: handleKeyDown,
        },
          h('span', { className: 'sc-step-dot sc-step-dot-' + status }, String(props.index)),
          h('span', { className: 'sc-step-title' }, T(step.titleKey)),
          h('span', { className: 'sc-step-status sc-status-' + status }, T(statusKey)),
          h('span', { className: 'sc-step-caret', 'aria-hidden': 'true' }, expanded ? '▾' : '▸')),
        body,
        h('div', { className: 'sc-step-foot' },
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-primary sc-btn-sm',
            disabled: busy || disabled,
            onClick: handleRun,
          }, busy ? T('running') : (status === 'done' ? T('rerun') : T('run'))),
          props.summaryText
            ? h('span', { className: 'sc-step-summary', title: props.summaryText }, props.summaryText)
            : null),
        localError || props.error
          ? h('div', { className: 'sc-inline-err' }, localError || props.error)
          : null)
    }

    function PipelinePanel(props) {
      var runningIds = props.runningStepIds || []
      var errors = props.errors || {}
      var cards = STEPS.map(function (step, index) {
        var status = stepStatusOf(step.id, props.doneSteps, runningIds)
        var summary = props.summaries && props.summaries[step.id] ? props.summaries[step.id] : null
        return h(StepCard, {
          key: (props.activeId || 'none') + ':' + step.id,
          index: index + 1,
          step: step,
          status: status,
          busy: runningIds.indexOf(step.id) >= 0,
          disabled: !props.activeId,
          summaryText: summaryTextOf(summary),
          error: errors['step:' + step.id] || '',
          annotationMapping: props.annotationMapping || {},
          annotationCount: props.annotationCount || 0,
          onRun: function (params) { props.onRunStep(step.id, params) },
        })
      })
      return h('aside', { className: 'sc-left', 'aria-label': T('pipeline') },
        h('div', { className: 'sc-panel-head' },
          h('span', { className: 'sc-panel-title' }, T('pipeline')),
          h('span', { className: 'sc-panel-note' },
            props.activeId ? (props.doneSteps.length + '/' + STEPS.length) : T('noData'))),
        h('div', { className: 'sc-panel-hint' }, T('pipelineHint')),
        h('div', { className: 'sc-steps' }, cards))
    }

    // ==========================================================================
    // 右栏：Agent 面板
    // ==========================================================================

    function AgentPanel(props) {
      var engine = props.engine || {}
      var badge = engineBadge(engine)
      var dataset = props.activeDataset
      var suggestions = props.suggestions || EMPTY_LIST
      var activity = props.activity || EMPTY_LIST
      var doneSteps = props.doneSteps || EMPTY_LIST
      var agentSuggestions = props.agentSuggestions || EMPTY_LIST
      var agentResults = props.agentResults || {}
      var agentInput = props.agentInput || ''

      var suggestionNodes = suggestions.map(function (suggestion) {
        var copyKey = 'sug:' + suggestion.id
        return h('div', { key: copyKey, className: 'sc-sug' },
          suggestion.title ? h('div', { className: 'sc-sug-title' }, suggestion.title) : null,
          suggestion.detail ? h('div', { className: 'sc-sug-detail' }, suggestion.detail) : null,
          h('div', { className: 'sc-sug-actions' },
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              onClick: function () { props.onCopy(copyKey, buildAgentInstruction(suggestion)) },
            }, props.copiedKey === copyKey ? T('copied') : T('copyInstruction'))))
      })

      var parsedNodes = agentSuggestions.map(function (suggestion) {
        var busy = props.agentBusyId === suggestion.id
        var result = agentResults[suggestion.id] || ''
        return h('div', { key: 'agent-sug-' + suggestion.id, className: 'sc-agent-sug' },
          h('div', { className: 'sc-sug-title' }, suggestion.title),
          suggestion.detail ? h('div', { className: 'sc-sug-detail' }, suggestion.detail) : null,
          h('div', { className: 'sc-agent-sug-actions' },
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm sc-btn-primary',
              disabled: Boolean(props.agentBusyId),
              onClick: function () { props.onAgentRun(suggestion) },
            }, busy ? T('agentRunning') : T('agentRun')),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: Boolean(props.agentBusyId),
              onClick: function () { props.onAgentCopy(suggestion) },
            }, T('copyInstruction'))),
          result
            ? h('div', { className: 'sc-agent-result' }, T('agentResult') + '：' + result)
            : null)
      })

      var activityNodes = activity.slice(0, 10).map(function (entry, index) {
        var meta = formatClock(entry.at) + (entry.ms !== null && entry.ms !== undefined ? ' · ' + formatDuration(entry.ms) : '')
        return h('div', { key: 'act-' + index, className: 'sc-act' },
          h('span', { className: 'sc-act-icon', 'aria-hidden': 'true' }, activityIcon(entry.actor)),
          h('div', { className: 'sc-act-body' },
            h('div', { className: 'sc-act-line' },
              h('span', { className: 'sc-act-action' }, entry.action),
              entry.ok ? null : h('span', { className: 'sc-bad' }, T('resultFail'))),
            entry.detail ? h('div', { className: 'sc-act-detail' }, entry.detail) : null,
            h('div', { className: 'sc-act-meta' }, meta)))
      })

      var askNode = h('div', { className: 'sc-agent-ask' },
        h('div', { className: 'sc-section-title' }, h('span', null, T('agentSend'))),
        h('textarea', {
          value: agentInput,
          placeholder: T('agentAskPlaceholder'),
          rows: 2,
          onChange: function (event) { props.onAgentInput(event.target.value) },
          onKeyDown: function (event) {
            if (event.key === 'Enter' && !event.shiftKey) {
              event.preventDefault()
              props.onAgentSubmit()
            }
          },
        }),
        h('div', { className: 'sc-agent-ask-actions' },
          h('span', { className: 'sc-muted' },
            agentSuggestions.length ? tf(T, 'agentParsed', agentSuggestions.length) : ''),
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-sm sc-btn-primary',
            disabled: !String(agentInput || '').trim(),
            onClick: props.onAgentSubmit,
          }, T('agentSend'))))

      return h('aside', { className: 'sc-right', 'aria-label': T('agentPanel') },
        h('div', { className: 'sc-panel-head' },
          h('span', { className: 'sc-panel-title' }, T('agentPanel'))),
        h('div', { className: 'sc-agent' },
          h('div', { className: 'sc-agent-status' },
            h('div', { className: 'sc-agent-row' },
              h('span', { className: 'sc-agent-k' }, T('statusTitle')),
              h('span', { className: 'sc-badge sc-badge-' + badge.tone }, T(badge.labelKey))),
            dataset
              ? h('div', { className: 'sc-agent-name' }, dataset.name || dataset.id)
              : h('div', { className: 'sc-muted' }, T('noData')),
            dataset
              ? h('div', { className: 'sc-agent-meta' },
                  formatNumber(dataset.nObs) + ' ' + T('cellsUnit') + ' × ' + formatNumber(dataset.nVars) + ' ' + T('genesUnit'))
              : null,
            h('div', { className: 'sc-agent-meta' }, T('engine') + ' · ' + T(badge.labelKey) + (engine.python ? ' · Python ' + engine.python : '')),
            dataset
              ? h('div', { className: 'sc-agent-meta' }, T('progress') + ' · ' + doneSteps.length + '/' + STEPS.length)
              : null,
            engine.error ? h('div', { className: 'sc-inline-err' }, engine.error) : null,
            (dataset && dataset.warnings && dataset.warnings.length)
              ? h('div', { className: 'sc-warn-text' }, dataset.warnings[0])
              : null),
          askNode,
          props.agentFallback
            ? h('div', { className: 'sc-sug' },
                h('div', { className: 'sc-sug-title' }, T('agentFallbackTitle')),
                h('div', { className: 'sc-sug-detail' }, T('agentFallbackHint')))
            : null,
          parsedNodes.length
            ? h('div', { className: 'sc-suggest' }, parsedNodes)
            : null,
          h('div', { className: 'sc-section-title' },
            h('span', null, T('suggestions')),
            h('span', { className: 'sc-badge' }, String(suggestions.length))),
          suggestions.length ? h('div', { className: 'sc-suggest' }, suggestionNodes) : Notice({ text: T('noSuggestions') }),
          props.copyError ? h('div', { className: 'sc-inline-err' }, props.copyError) : null,
          h('div', { className: 'sc-section-title' }, h('span', null, T('activity'))),
          activity.length ? h('div', { className: 'sc-timeline' }, activityNodes) : Notice({ text: T('activityEmpty') }),
          h('div', { className: 'sc-agent-hint' }, T('agentHint'))))
    }

    // ==========================================================================
    // 中栏 tabs：嵌入图 / QC / 标记基因 / 聚类 / 记录 / 日志
    // ==========================================================================

    var EMPTY_LIST = []

    var TAB_DEFS = [
      { id: 'embedding', labelKey: 'tabEmbedding' },
      { id: 'qc', labelKey: 'tabQc' },
      { id: 'markers', labelKey: 'tabMarkers' },
      { id: 'clusters', labelKey: 'tabClusters' },
      { id: 'record', labelKey: 'tabRecord' },
      { id: 'log', labelKey: 'tabLog' },
      { id: 'analysis', labelKey: 'tabAnalysis' },
    ]

    function Fragment(props) {
      return props.children
    }

    function TabsBar(props) {
      return h('div', { className: 'sc-tabs', role: 'tablist' }, TAB_DEFS.map(function (tab) {
        return h('button', {
          key: tab.id,
          type: 'button',
          role: 'tab',
          'aria-selected': props.tab === tab.id ? 'true' : 'false',
          className: 'sc-tab' + (props.tab === tab.id ? ' sc-tab-on' : ''),
          onClick: function () { props.onChange(tab.id) },
        }, T(tab.labelKey))
      }))
    }

    function EmptyState(props) {
      var busy = props.busy || {}
      return h('div', { className: 'sc-empty' },
        h('div', { className: 'sc-empty-card' },
          h('div', { className: 'sc-empty-title' }, T('emptyTitle')),
          h('div', { className: 'sc-empty-hint' }, T('emptyHint')),
          h('div', { className: 'sc-empty-actions' },
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-primary',
              disabled: Boolean(busy.import),
              onClick: props.onOpenImport,
            }, T('importData')),
            h('button', {
              type: 'button',
              className: 'sc-btn',
              disabled: Boolean(busy.demo),
              onClick: props.onDemo,
            }, busy.demo ? T('running') : T('loadDemo')))))
    }

    // ==========================================================================
    // 嵌入图：canvas 交互散点（滚轮缩放以鼠标为锚 / 拖拽平移 / 双击复位 / 悬停 tooltip）
    // ==========================================================================

    function EmbeddingCanvas(props) {
      var canvasRef = React.useRef(null)
      var wrapRef = React.useRef(null)
      var viewRef = React.useRef(null)
      var sizeRef = React.useRef({ width: 0, height: 0 })
      var dragRef = React.useRef(null)
      var tickState = React.useState(0)
      var tick = tickState[0]
      var setTick = tickState[1]
      var sizeState = React.useState({ width: 0, height: 0 })
      var size = sizeState[0]
      var setSize = sizeState[1]
      var tooltipState = React.useState(null)
      var tooltip = tooltipState[0]
      var setTooltip = tooltipState[1]
      var data = props.data

      // 尺寸测量（ResizeObserver 优先）
      React.useEffect(function () {
        var element = wrapRef.current
        if (!element || typeof element.getBoundingClientRect !== 'function') return undefined
        function measure() {
          var rect = element.getBoundingClientRect()
          var width = Math.max(1, Math.round(rect.width))
          var height = Math.max(1, Math.round(rect.height))
          sizeRef.current = { width: width, height: height }
          setSize(function (previous) {
            return previous.width === width && previous.height === height ? previous : { width: width, height: height }
          })
        }
        measure()
        var observer = null
        if (typeof ResizeObserver === 'function') {
          observer = new ResizeObserver(measure)
          observer.observe(element)
        } else if (typeof window !== 'undefined' && window.addEventListener) {
          window.addEventListener('resize', measure)
        }
        return function () {
          if (observer) {
            observer.disconnect()
          } else if (typeof window !== 'undefined' && window.removeEventListener) {
            window.removeEventListener('resize', measure)
          }
        }
      }, [])

      // 数据换了：重新适配视图（缩放/平移只在交互里改，不重新请求）
      React.useEffect(function () {
        if (!data || !data.x.length) {
          viewRef.current = null
          return
        }
        var width = sizeRef.current.width || 640
        var height = sizeRef.current.height || 380
        viewRef.current = fitView(computeBounds(data.x, data.y), width, height, 26)
        setTick(function (value) { return value + 1 })
      }, [data])

      // 首次拿到尺寸时补一次适配
      React.useEffect(function () {
        if (!data || !data.x.length || !size.width || !size.height || viewRef.current) return
        viewRef.current = fitView(computeBounds(data.x, data.y), size.width, size.height, 26)
        setTick(function (value) { return value + 1 })
      }, [data, size])

      // 原生 wheel：React 的 onWheel 在宿主里可能是 passive，preventDefault 会失效
      React.useEffect(function () {
        var canvas = canvasRef.current
        if (!canvas || typeof canvas.addEventListener !== 'function') return undefined
        function handleWheel(event) {
          var view = viewRef.current
          if (!view) return
          if (typeof event.preventDefault === 'function') event.preventDefault()
          var rect = canvas.getBoundingClientRect()
          var local = { x: event.clientX - rect.left, y: event.clientY - rect.top }
          var factor = event.deltaY < 0 ? 1.15 : 1 / 1.15
          var currentSize = sizeRef.current
          viewRef.current = zoomViewAt(
            view,
            local,
            factor,
            { x: currentSize.width / 2, y: currentSize.height / 2 },
            0.02,
            500,
          )
          setTooltip(null)
          setTick(function (value) { return value + 1 })
        }
        canvas.addEventListener('wheel', handleWheel, { passive: false })
        return function () { canvas.removeEventListener('wheel', handleWheel) }
      }, [])

      function resetView() {
        if (!data || !data.x.length) return
        var currentSize = sizeRef.current
        if (!currentSize.width || !currentSize.height) return
        viewRef.current = fitView(computeBounds(data.x, data.y), currentSize.width, currentSize.height, 26)
        setTooltip(null)
        setTick(function (value) { return value + 1 })
      }

      function localPoint(event) {
        var canvas = canvasRef.current
        if (!canvas || !event || typeof canvas.getBoundingClientRect !== 'function') return null
        var rect = canvas.getBoundingClientRect()
        return { x: event.clientX - rect.left, y: event.clientY - rect.top }
      }

      function handleMouseDown(event) {
        if (event.button !== 0) return
        var local = localPoint(event)
        if (!local) return
        dragRef.current = { x: local.x, y: local.y }
        if (typeof event.preventDefault === 'function') event.preventDefault()
      }

      function handleMouseUp() {
        dragRef.current = null
      }

      function handleMouseLeave() {
        dragRef.current = null
        setTooltip(null)
      }

      function handleMouseMove(event) {
        var local = localPoint(event)
        var view = viewRef.current
        if (!local || !view) return
        var drag = dragRef.current
        if (drag) {
          var deltaX = local.x - drag.x
          var deltaY = local.y - drag.y
          drag.x = local.x
          drag.y = local.y
          viewRef.current = panView(view, deltaX, deltaY)
          setTooltip(null)
          setTick(function (value) { return value + 1 })
          return
        }
        var currentSize = sizeRef.current
        var center = { x: currentSize.width / 2, y: currentSize.height / 2 }
        var index = findNearestPoint(data.x, data.y, view, center, local, 14)
        if (index < 0) {
          if (tooltip) setTooltip(null)
          return
        }
        var label = ''
        if (!data.continuous) {
          var group = data.groupIndex[index]
          var category = data.categories[group]
          label = category ? category.label : (group === undefined ? '' : String(group))
        }
        var text = tooltipText({
          index: index,
          continuous: data.continuous,
          title: data.legendTitle,
          value: data.continuous ? data.values[index] : undefined,
          label: label,
        })
        var tipX = local.x + 12
        var tipY = local.y - 30
        if (tipX > currentSize.width - 170) tipX = Math.max(4, local.x - 160)
        if (tipY < 4) tipY = local.y + 16
        setTooltip({ x: tipX, y: tipY, text: text })
      }

      // 绘制：处理 devicePixelRatio；缩放/平移只改坐标变换，不重新请求数据。
      React.useEffect(function () {
        var canvas = canvasRef.current
        if (!canvas || !data || !data.x.length || !size.width || !size.height) return
        var dpr = (typeof window !== 'undefined' && window.devicePixelRatio) || 1
        canvas.width = Math.max(1, Math.round(size.width * dpr))
        canvas.height = Math.max(1, Math.round(size.height * dpr))
        var context = null
        try {
          context = canvas.getContext('2d')
        } catch (error) {
          context = null
        }
        if (!context) return
        context.setTransform(dpr, 0, 0, dpr, 0, 0)
        context.clearRect(0, 0, size.width, size.height)
        var view = viewRef.current
        if (!view) {
          view = fitView(computeBounds(data.x, data.y), size.width, size.height, 26)
          viewRef.current = view
        }
        var center = { x: size.width / 2, y: size.height / 2 }
        var pointSize = clamp(Math.sqrt((size.width * size.height) / Math.max(1, data.x.length)) * 0.55, 2.2, 5.5)
        var half = pointSize / 2
        var colors = new Array(data.x.length)
        var index
        if (data.continuous) {
          var min = Infinity
          var max = -Infinity
          for (index = 0; index < data.values.length; index++) {
            var value = data.values[index]
            if (value < min) min = value
            if (value > max) max = value
          }
          for (index = 0; index < data.x.length; index++) {
            colors[index] = pickContinuousColor(data.values[index], min, max)
          }
        } else {
          var groupColors = []
          for (index = 0; index < data.categories.length; index++) {
            groupColors[index] = data.categories[index].color || PALETTE[index % PALETTE.length]
          }
          for (index = 0; index < data.x.length; index++) {
            colors[index] = groupColors[data.groupIndex[index]] || PALETTE[0]
          }
        }
        var highlight = props.highlight
        var dimming = highlight !== null && highlight !== undefined && !data.continuous
        for (index = 0; index < data.x.length; index++) {
          var projected = worldToScreen({ x: data.x[index], y: data.y[index] }, view, center)
          if (projected.x < -pointSize || projected.y < -pointSize || projected.x > size.width + pointSize || projected.y > size.height + pointSize) continue
          if (dimming) {
            context.globalAlpha = String(data.groupIndex[index]) === String(highlight) ? 1 : 0.15
          } else {
            context.globalAlpha = 1
          }
          context.fillStyle = colors[index]
          context.fillRect(projected.x - half, projected.y - half, pointSize, pointSize)
        }
        context.globalAlpha = 1
      }, [data, props.highlight, size, tick])

      var legendNode = null
      if (data) {
        if (data.continuous) {
          var minValue = Infinity
          var maxValue = -Infinity
          for (var valueIndex = 0; valueIndex < data.values.length; valueIndex++) {
            var item = data.values[valueIndex]
            if (item < minValue) minValue = item
            if (item > maxValue) maxValue = item
          }
          if (!Number.isFinite(minValue)) minValue = 0
          if (!Number.isFinite(maxValue)) maxValue = 0
          legendNode = h('div', { className: 'sc-gradient' },
            h('span', null, formatNumber(minValue, 2)),
            h('span', { className: 'sc-gradient-bar', style: { background: 'linear-gradient(90deg,' + CONTINUOUS_STOPS.join(',') + ')' } }),
            h('span', null, formatNumber(maxValue, 2)),
            h('span', { className: 'sc-muted' }, data.legendTitle || ''))
        } else if (data.categories.length) {
          legendNode = h('div', { className: 'sc-legend' },
            h('div', { className: 'sc-muted' }, T('legend') + (data.legendTitle ? ' · ' + data.legendTitle : '')),
            h('div', { className: 'sc-legend-row' },
              data.categories.map(function (category) {
                var on = props.highlight !== null && props.highlight !== undefined && String(props.highlight) === category.id
                return h('button', {
                  key: category.id,
                  type: 'button',
                  className: 'sc-legend-chip' + (on ? ' sc-legend-chip-on' : ''),
                  title: category.label,
                  onClick: function () { props.onHighlight(on ? null : category.id) },
                },
                  h('span', { className: 'sc-legend-dot', style: { background: category.color } }),
                  h('span', null, category.label),
                  category.n ? h('span', { className: 'sc-legend-n' }, formatNumber(category.n)) : null)
              })))
        }
      }

      var metaText = (data.layout ? data.layout + ' · ' : '') + formatNumber(data.nPoints) + ' ' + T('pointsUnit') +
        (data.downsampled ? ' · ' + tf(T, 'downsampleNotice', formatNumber(data.nPoints)) : '')

      return h('div', { className: 'sc-canvas-stack' },
        h('div', { className: 'sc-canvas-wrap', ref: wrapRef },
          h('canvas', {
            ref: canvasRef,
            className: 'sc-canvas',
            onMouseDown: handleMouseDown,
            onMouseMove: handleMouseMove,
            onMouseUp: handleMouseUp,
            onMouseLeave: handleMouseLeave,
            onDoubleClick: resetView,
          }),
          tooltip
            ? h('div', { className: 'sc-tip', style: { left: tooltip.x + 'px', top: tooltip.y + 'px' } }, tooltip.text)
            : null,
          h('button', { type: 'button', className: 'sc-btn sc-btn-sm sc-canvas-reset', onClick: resetView }, T('resetView')),
          h('div', { className: 'sc-canvas-meta' }, metaText + ' · ' + T('hoverHint'))),
        legendNode)
    }

    function colorValueOf(colorBy) {
      if (!colorBy || typeof colorBy !== 'object') return 'cluster:clusters'
      if (colorBy.kind === 'gene') return 'gene'
      return String(colorBy.kind || 'cluster') + ':' + String(colorBy.key || 'clusters')
    }

    function EmbeddingTab(props) {
      var colorValue = colorValueOf(props.colorBy)
      var geneMode = colorValue === 'gene'
      var genes = props.genes || []
      var colorOptions = props.colorOptions && props.colorOptions.length
        ? props.colorOptions
        : buildColorOptions(EMPTY_LIST)
      var geneDatalist = h('datalist', { id: 'sc-gene-list' }, genes.map(function (gene) {
        return h('option', { key: gene, value: gene })
      }))
      var tools = h('div', { className: 'sc-embed-tools' },
        h('label', { className: 'sc-tool-field' },
          h('span', null, T('colorBy')),
          h('select', {
            className: 'sc-input sc-select',
            value: colorValue,
            disabled: props.loading,
            onChange: function (event) { props.onColorChange(event.target.value) },
          }, colorOptions.map(function (option) {
            return h('option', { key: option.value, value: option.value },
              option.labelKey ? T(option.labelKey) : (option.label || option.value))
          }))),
        geneMode
          ? h('div', { className: 'sc-gene-form' },
              h('input', {
                className: 'sc-input sc-gene-input',
                type: 'text',
                list: 'sc-gene-list',
                value: props.geneQuery,
                placeholder: T('genePlaceholder'),
                onChange: function (event) { props.onGeneQuery(event.target.value) },
                onKeyDown: function (event) {
                  if (event.key === 'Enter') {
                    event.preventDefault()
                    props.onGeneSubmit()
                  }
                },
              }),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm',
                disabled: props.loading || !String(props.geneQuery || '').trim(),
                onClick: props.onGeneSubmit,
              }, T('geneApply')),
              geneDatalist)
          : geneDatalist,
        h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: props.loading, onClick: props.onRefresh }, T('refresh')),
        props.layout ? h('span', { className: 'sc-muted' }, props.layout) : null)

      var body = null
      if (!props.hasEmbed) {
        body = h('div', { className: 'sc-guide' },
          h('div', { className: 'sc-guide-title' }, T('noEmbeddingTitle')),
          h('div', { className: 'sc-guide-hint' }, T('noEmbeddingHint')),
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-primary',
            disabled: !props.pcaDone || props.running,
            onClick: props.onRunEmbed,
          }, props.running ? T('running') : T('runEmbed')),
          !props.pcaDone ? h('div', { className: 'sc-hint-line' }, T('embedNeedPca')) : null)
      } else if (geneMode && !String(props.colorBy.key || '').trim()) {
        body = h('div', { className: 'sc-guide' },
          h('div', { className: 'sc-guide-hint' }, T('geneEmptyHint')))
      } else if (props.loading && !props.data) {
        body = LoadingLine({ text: T('loadingEmbedding') })
      } else if (props.error) {
        body = ErrorBanner({ message: props.error, onRetry: props.onRefresh, retryLabel: T('retry') })
      } else if (!props.data) {
        body = Notice({ text: T('embedEmpty') })
      } else {
        body = h(EmbeddingCanvas, {
          data: props.data,
          highlight: props.highlight,
          onHighlight: props.onHighlight,
        })
      }

      return h(Fragment, null, tools, body)
    }

    // ==========================================================================
    // QC tab
    // ==========================================================================

    var QC_METRIC_DEFS = [
      { id: 'n_genes', label: 'n_genes' },
      { id: 'n_counts', label: 'n_counts' },
      { id: 'pct_counts_mt', label: 'pct_counts_mt' },
    ]

    function QcTab(props) {
      if (!props.activeId) return Notice({ text: T('noData') })
      var metrics = props.metrics || []
      return h(Fragment, null,
        h('div', { className: 'sc-toolbar' },
          h('span', { className: 'sc-muted' }, T('qcMetrics')),
          QC_METRIC_DEFS.map(function (metric) {
            return h('label', { key: metric.id, className: 'sc-check' },
              h('input', {
                type: 'checkbox',
                checked: metrics.indexOf(metric.id) >= 0,
                disabled: props.loading,
                onChange: function () { props.onToggleMetric(metric.id) },
              }),
              metric.label)
          }),
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-primary',
            disabled: props.loading || !metrics.length,
            onClick: props.onGenerate,
          }, props.loading ? T('running') : T('generateQc'))),
        props.error
          ? ErrorBanner({ message: props.error, onRetry: props.onGenerate, retryLabel: T('retry') })
          : null,
        props.figure
          ? h('div', { className: 'sc-figure' },
              h('img', { src: props.figure.url, alt: T('qcFigureAlt'), loading: 'lazy' }),
              h('div', { className: 'sc-figure-meta' }, (props.figure.width || 0) + '×' + (props.figure.height || 0) + 'px'))
          : h('div', { className: 'sc-guide' },
              h('div', { className: 'sc-guide-hint' }, T('qcEmptyHint'))),
        !props.qcDone ? h('div', { className: 'sc-hint-line' }, T('qcNeedStep')) : null)
    }

    // ==========================================================================
    // 标记基因 tab
    // ==========================================================================

    // ==========================================================================
    // v6：降维聚类控制卡（聚类 tab 顶部）
    // ==========================================================================

    function DimredCard(props) {
      var dim = props.dim || DIMRED_DEFAULTS
      var batchOptions = props.batchOptions || EMPTY_LIST
      var busy = Boolean(props.runningRecipe || props.runningRecluster)
      var warnings = props.warnings || EMPTY_LIST
      return h('section', { className: 'sc-dimred' },
        h('div', { className: 'sc-dimred-head' },
          h('div', null,
            h('div', { className: 'sc-dimred-title' }, T('dimredTitle')),
            h('div', { className: 'sc-dimred-sub' }, T('dimredSub'))),
          h('div', { className: 'sc-dimred-actions' },
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-primary',
              disabled: busy || !props.activeId,
              onClick: props.onRunRecipe,
            }, props.runningRecipe ? T('dimredRunBusy') : T('dimredRun')),
            h('button', {
              type: 'button',
              className: 'sc-btn',
              disabled: busy || !props.activeId,
              onClick: props.onRecluster,
            }, props.runningRecluster ? T('dimredReclusterBusy') : T('dimredRecluster')))),
        h('div', { className: 'sc-dimred-grid' },
          h('label', { className: 'sc-dimred-field' },
            h('span', null, T('dimredEngine')),
            h('select', {
              className: 'sc-input sc-select',
              value: dim.engine,
              disabled: busy,
              onChange: function (event) { props.onEngine(event.target.value) },
            },
              h('option', { value: 'seurat' }, T('dimredEngineSeurat')),
              h('option', { value: 'builtin' }, T('dimredEngineBuiltin')))),
          h('label', { className: 'sc-check sc-dimred-field' },
            h('input', {
              type: 'checkbox',
              checked: dim.harmony === true,
              disabled: busy || dim.engine !== 'seurat',
              onChange: function (event) { props.onHarmony(event.target.checked) },
            }),
            T('dimredHarmony')),
          h('label', { className: 'sc-dimred-field' },
            h('span', null, T('dimredBatchKey')),
            h('select', {
              className: 'sc-input sc-select',
              value: dim.batchKey || '',
              disabled: busy || dim.engine !== 'seurat',
              title: T('dimredBatchHint'),
              onChange: function (event) { props.onBatchKey(event.target.value) },
            },
              h('option', { value: '' }, T('dimredBatchAuto')),
              batchOptions.map(function (column) {
                return h('option', { key: 'batch-' + column, value: column }, column)
              }))),
          h('label', { className: 'sc-dimred-field' },
            h('span', null, T('dimredNPcs')),
            h('input', {
              className: 'sc-input',
              type: 'number',
              min: 2,
              max: 200,
              step: 1,
              value: String(dim.nPcs),
              disabled: busy,
              onChange: function (event) { props.onNPcs(event.target.value) },
            })),
          h('label', { className: 'sc-dimred-field' },
            h('span', null, T('dimredNeighbors')),
            h('input', {
              className: 'sc-input',
              type: 'number',
              min: 2,
              max: 200,
              step: 1,
              value: String(dim.neighbors),
              disabled: busy,
              onChange: function (event) { props.onNeighbors(event.target.value) },
            })),
          h('label', { className: 'sc-dimred-field' },
            h('span', null, T('dimredTopGenes')),
            h('input', {
              className: 'sc-input',
              type: 'number',
              min: 100,
              max: 10000,
              step: 100,
              value: String(dim.nTopGenes),
              disabled: busy,
              onChange: function (event) { props.onTopGenes(event.target.value) },
            })),
          h('div', { className: 'sc-dimred-slider-row' },
            h('span', null, T('dimredResolution')),
            h('input', {
              className: 'sc-dimred-slider',
              type: 'range',
              min: 0.05,
              max: 1.5,
              step: 0.05,
              value: String(dim.resolution),
              disabled: busy,
              onChange: function (event) { props.onResolution(event.target.value) },
            }),
            h('span', { className: 'sc-dimred-value' }, String(dim.resolution)),
            h('span', { className: 'sc-muted' }, T('dimredHint')))),
        dim.engine === 'seurat'
          ? h('div', { className: 'sc-hint-line' }, T('dimredSeuratHint'))
          : null,
        props.error ? h('div', { className: 'sc-inline-err' }, props.error) : null,
        warnings.length
          ? h('div', { className: 'sc-dimred-warns' },
              h('div', { className: 'sc-muted' }, T('dimredWarnings')),
              warnings.map(function (warning, index) {
                return h('div', { key: 'dw-' + index }, warning)
              }))
          : null)
    }

    // ==========================================================================
    // v6：手动注释工作台（标记基因 tab 顶部）
    // ==========================================================================

    function AnnotationWorkbench(props) {
      var groups = props.groups || EMPTY_LIST
      var annotation = props.annotation || {}
      var rows = buildClusterAnnotationRows(groups, annotation)
      var selectedCluster = props.selectedCluster === null || props.selectedCluster === undefined ? '' : String(props.selectedCluster)
      var selectedGene = String(props.selectedGene || '')
      var selectedGroup = null
      for (var index = 0; index < groups.length; index++) {
        if (String(groups[index].id) === selectedCluster) selectedGroup = groups[index]
      }
      var genes = selectedGroup && Array.isArray(selectedGroup.genes) ? selectedGroup.genes : []
      var clusterButtons = rows.map(function (row) {
        var on = row.id === selectedCluster
        return h('button', {
          key: 'annot-cluster-' + row.id,
          type: 'button',
          className: 'sc-annot-cluster' + (on ? ' sc-annot-cluster-on' : ''),
          onClick: function () { props.onSelectCluster(row.id) },
        },
          h('span', null, T('annotClusters') + ' ' + row.label + (row.nCells ? ' · ' + formatNumber(row.nCells) : '')),
          h('span', { className: 'sc-annot-cluster-top' }, row.topGenes.length ? row.topGenes.join(' · ') : '—'))
      })
      var geneChips = genes.map(function (gene, geneIndex) {
        var name = gene && gene.gene ? String(gene.gene) : ('gene-' + geneIndex)
        return h('button', {
          key: 'annot-gene-' + name,
          type: 'button',
          className: 'sc-annot-gene-chip' + (name === selectedGene ? ' sc-annot-gene-chip-on' : ''),
          onClick: function () { props.onSelectGene(name) },
        }, name)
      })
      var markerRows = genes.map(function (gene, geneIndex) {
        var name = gene && gene.gene ? String(gene.gene) : ''
        return h('tr', { key: 'annot-row-' + (name || geneIndex) },
          h('td', null, h('button', {
            type: 'button',
            className: 'sc-annot-gene-chip' + (name === selectedGene ? ' sc-annot-gene-chip-on' : ''),
            onClick: function () { props.onSelectGene(name) },
          }, name || '—')),
          h('td', { className: 'sc-td-num' }, formatNumber(gene.logfc, 2)),
          h('td', { className: 'sc-td-num' }, formatPValue(gene.pvalAdj)),
          h('td', { className: 'sc-td-num' }, formatPct(gene.pctIn)))
      })
      var mappingChips = rows.filter(function (row) { return row.cellType }).map(function (row) {
        return h('span', { key: 'annot-map-' + row.id, className: 'sc-annot-map-chip' }, 'cluster ' + row.label + ' → ' + row.cellType)
      })
      var clusterInput = selectedCluster !== '' && selectedCluster !== null
        ? h('div', { className: 'sc-annot-actions' },
            h('span', { className: 'sc-muted' }, T('annotClusters') + ' ' + selectedCluster + ' →'),
            h('input', {
              className: 'sc-input sc-sample-name',
              type: 'text',
              value: annotation[selectedCluster] || '',
              placeholder: T('cellTypePlaceholder'),
              disabled: Boolean(props.saving),
              onChange: function (event) { props.onChangeAnnotation(selectedCluster, event.target.value) },
            }),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-primary sc-btn-sm',
              disabled: Boolean(props.saving) || !Object.keys(annotation).length,
              onClick: props.onSaveAnnotation,
            }, props.saving ? T('running') : T('annotSave')),
            props.savedMessage ? h('span', { className: 'sc-ok' }, props.savedMessage) : null)
        : h('div', { className: 'sc-hint-line' }, T('samplesGroupAdderHint'))

      return h('section', { className: 'sc-annot' },
        h('div', { className: 'sc-annot-head' },
          h('div', null,
            h('div', { className: 'sc-annot-title' }, T('annotTitle')),
            h('div', { className: 'sc-record-note' }, T('annotSubtitle'))),
          h('div', { className: 'sc-annot-actions' },
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: Boolean(props.csvOpen),
              onClick: props.onToggleCsv,
            }, T('annotCsv')))),
        h('div', { className: 'sc-annot-layout' },
          h('div', { className: 'sc-annot-clusters' },
            rows.length ? clusterButtons : h('div', { className: 'sc-record-note' }, T('annotNeedMarkers'))),
          h('div', { className: 'sc-annot-main' },
            h('div', { className: 'sc-inline-row' },
              h('span', { className: 'sc-muted' }, T('annotSelectedGene') + '：' + (selectedGene || '—')),
              h('span', { className: 'sc-muted' }, T('annotNoGeneHint'))),
            geneChips.length ? h('div', { className: 'sc-annot-genes' }, geneChips) : null,
            h('div', { className: 'sc-annot-actions' },
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm',
                disabled: Boolean(props.plotBusy) || !selectedGene,
                onClick: function () { props.onPlot('dotplot') },
              }, props.plotBusy ? T('annotPlotBusy') : T('annotDotplot')),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm',
                disabled: Boolean(props.plotBusy) || !selectedGene,
                onClick: function () { props.onPlot('heatmap') },
              }, T('annotHeatmap')),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm',
                disabled: Boolean(props.plotBusy) || !selectedGene,
                onClick: function () { props.onFeaturePlot() },
              }, T('annotFeatureplot'))),
            h('div', { className: 'sc-hint-line' }, T('annotPlotNote')),
            markerRows.length
              ? h('div', { className: 'sc-table-wrap' },
                  h('table', { className: 'sc-table' },
                    h('thead', null, h('tr', null,
                      h('th', null, T('annotGene')),
                      h('th', { className: 'sc-td-num' }, T('tableLogfc')),
                      h('th', { className: 'sc-td-num' }, T('tablePval')),
                      h('th', { className: 'sc-td-num' }, T('tablePctIn')))),
                    h('tbody', null, markerRows)))
              : null,
            clusterInput,
            props.plotError ? h('div', { className: 'sc-inline-err' }, props.plotError) : null,
            props.plotFigure
              ? h('div', { className: 'sc-figure sc-figure-sm' },
                  h('img', { src: props.plotFigure.url, alt: T('annotMarkerTable'), loading: 'lazy' }))
              : null,
            mappingChips.length
              ? h('div', { className: 'sc-annot-mapping' },
                  h('div', { className: 'sc-muted' }, T('annotMappingTable')),
                  mappingChips)
              : null,
            props.saveError ? h('div', { className: 'sc-inline-err' }, props.saveError) : null,
            props.csvOpen
              ? h('div', { className: 'sc-annot-csv' },
                  h('textarea', {
                    value: props.csvText || '',
                    placeholder: T('annotCsvPlaceholder'),
                    onChange: function (event) { props.onCsvText(event.target.value) },
                  }),
                  h('div', { className: 'sc-annot-actions' },
                    h('button', {
                      type: 'button',
                      className: 'sc-btn sc-btn-sm sc-btn-primary',
                      disabled: !String(props.csvText || '').trim(),
                      onClick: props.onApplyCsv,
                    }, T('annotCsvApply')),
                    h('button', { type: 'button', className: 'sc-btn sc-btn-sm', onClick: props.onToggleCsv }, T('cancel'))))
              : null)))
    }

    function MarkersTab(props) {
      if (!props.activeId) return Notice({ text: T('noData') })
      if (!props.clusterDone) {
        return h('div', { className: 'sc-guide' },
          h('div', { className: 'sc-guide-hint' }, T('markersNeedCluster')))
      }
      var groups = props.groups || []
      var rows = []
      groups.forEach(function (group) {
        var genes = group.genes || []
        if (!genes.length) {
          rows.push(h('tr', { key: 'g-' + group.id },
            h('td', { className: 'sc-td-group' }, group.label),
            h('td', { colSpan: 4, className: 'sc-muted' }, '—')))
          return
        }
        genes.forEach(function (gene, index) {
          rows.push(h('tr', { key: group.id + '-' + (gene.gene || index) },
            index === 0
              ? h('td', { className: 'sc-td-group', rowSpan: genes.length },
                  group.label + (group.nCells ? ' (' + formatNumber(group.nCells) + ')' : ''))
              : null,
            h('td', null, gene.gene),
            h('td', { className: 'sc-td-num' }, formatNumber(gene.logfc, 2)),
            h('td', { className: 'sc-td-num' }, formatPValue(gene.pvalAdj)),
            h('td', { className: 'sc-td-num' }, formatPct(gene.pctIn))))
        })
      })
      var bodyNode = null
      if (props.loading && !groups.length) {
        bodyNode = LoadingLine({ text: T('loading') })
      } else if (props.error) {
        bodyNode = ErrorBanner({ message: props.error, onRetry: props.onReload, retryLabel: T('retry') })
      } else if (!groups.length) {
        bodyNode = Notice({ text: T('markersEmpty') })
      } else {
        bodyNode = h('div', { className: 'sc-table-wrap' },
          h('table', { className: 'sc-table' },
            h('thead', null, h('tr', null,
              h('th', null, T('tableCluster')),
              h('th', null, T('tableGene')),
              h('th', { className: 'sc-td-num' }, T('tableLogfc')),
              h('th', { className: 'sc-td-num' }, T('tablePval')),
              h('th', { className: 'sc-td-num' }, T('tablePctIn')))),
            h('tbody', null, rows)))
      }
      return h(Fragment, null,
        h(AnnotationWorkbench, {
          groups: groups,
          annotation: props.annotation,
          selectedCluster: props.selectedCluster,
          selectedGene: props.selectedGene,
          onSelectCluster: props.onSelectCluster,
          onSelectGene: props.onSelectGene,
          onChangeAnnotation: props.onChangeAnnotation,
          onSaveAnnotation: props.onSaveAnnotation,
          savedMessage: props.annotationSaved,
          saveError: props.annotationError,
          saving: props.savingAnnotation,
          plotBusy: props.plotBusy,
          plotError: props.plotError,
          plotFigure: props.plotFigure,
          onPlot: props.onPlotMarker,
          onFeaturePlot: props.onFeaturePlot,
          csvOpen: props.csvOpen,
          csvText: props.csvText,
          onToggleCsv: props.onToggleCsv,
          onCsvText: props.onCsvText,
          onApplyCsv: props.onApplyCsv,
        }),
        h('div', { className: 'sc-toolbar' },
          h('label', { className: 'sc-tool-field' },
            h('span', null, T('paramMarkerMethod')),
            h('select', {
              className: 'sc-input sc-select',
              value: props.method,
              disabled: props.loading,
              onChange: function (event) { props.onMethodChange(event.target.value) },
            },
              h('option', { value: 't-test' }, 't-test'),
              h('option', { value: 'wilcoxon' }, 'wilcoxon'))),
          h('label', { className: 'sc-tool-field' },
            h('span', null, T('paramTopN')),
            h('select', {
              className: 'sc-input sc-select',
              value: String(props.topN),
              disabled: props.loading,
              onChange: function (event) { props.onTopNChange(Number(event.target.value)) },
            },
              h('option', { value: '5' }, '5'),
              h('option', { value: '10' }, '10'),
              h('option', { value: '20' }, '20'))),
          h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: props.loading, onClick: props.onReload }, T('reloadMarkers')),
          h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: props.dotplotLoading, onClick: props.onDotplot },
            props.dotplotLoading ? T('running') : T('generateDotplot'))),
        bodyNode,
        props.dotplotError
          ? ErrorBanner({ message: props.dotplotError, onRetry: props.onDotplot, retryLabel: T('retry') })
          : null,
        props.dotplot
          ? h('div', { className: 'sc-figure sc-figure-sm' },
              h('img', { src: props.dotplot.url, alt: T('dotplotAlt'), loading: 'lazy' }))
          : null)
    }

    // ==========================================================================
    // 聚类 tab：cluster_bar + 细胞类型注释
    // ==========================================================================

    function ClustersTab(props) {
      if (!props.activeId) return Notice({ text: T('noData') })
      if (!props.clusterDone) {
        return h('div', { className: 'sc-guide' },
          h('div', { className: 'sc-guide-hint' }, T('markersNeedCluster')))
      }
      var groups = props.groups || []
      var listNode = null
      if (props.markersLoading && !groups.length) {
        listNode = LoadingLine({ text: T('loading') })
      } else if (props.markersError) {
        listNode = ErrorBanner({ message: props.markersError, onRetry: props.onReloadMarkers, retryLabel: T('retry') })
      } else if (!groups.length) {
        listNode = Notice({ text: T('markersEmpty') })
      } else {
        listNode = h('div', { className: 'sc-cluster-list' }, groups.map(function (group) {
          return h('div', { key: group.id, className: 'sc-cluster-row' },
            h('span', { className: 'sc-cluster-id' }, group.label),
            h('span', { className: 'sc-cluster-n' }, formatNumber(group.nCells) + ' ' + T('cellsUnit')),
            h('input', {
              className: 'sc-input sc-cluster-input',
              type: 'text',
              placeholder: T('cellTypePlaceholder'),
              value: props.annotation[group.id] || '',
              onChange: function (event) { props.onChangeAnnotation(group.id, event.target.value) },
            }))
        }))
      }
      return h(Fragment, null,
        props.dimredCard ? props.dimredCard : null,
        h('div', { className: 'sc-toolbar' },
          h('span', { className: 'sc-muted' }, T('clusterBarTitle')),
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-sm',
            disabled: props.barLoading,
            onClick: props.onGenerateBar,
          }, props.barLoading ? T('running') : T('generateClusterBar'))),
        props.barError
          ? ErrorBanner({ message: props.barError, onRetry: props.onGenerateBar, retryLabel: T('retry') })
          : null,
        props.barFigure
          ? h('div', { className: 'sc-figure sc-figure-sm' },
              h('img', { src: props.barFigure.url, alt: T('clusterBarTitle'), loading: 'lazy' }))
          : null,
        h('div', { className: 'sc-section-title' }, h('span', null, T('annotationTitle'))),
        h('div', { className: 'sc-hint-line' }, T('annotationHint')),
        listNode,
        h('div', { className: 'sc-toolbar' },
          h('button', {
            type: 'button',
            className: 'sc-btn sc-btn-primary',
            disabled: props.saving || !props.hasMapping,
            onClick: props.onSave,
          }, props.saving ? T('running') : T('saveAnnotation')),
          props.saved ? h('span', { className: 'sc-ok' }, props.saved) : null,
          props.saveError ? h('span', { className: 'sc-bad' }, props.saveError) : null))
    }

    // ==========================================================================
    // v4：样本 & 分组面板（记录 tab 顶部）
    // ==========================================================================

    function SamplesPanel(props) {
      var rows = props.rows || EMPTY_LIST
      var saving = Boolean(props.saving)
      var groupOptionsList = props.groupOptions || [{ value: '', label: '（未分组）' }]
      var autoPreview = props.autoPreview
      var importPreview = props.importPreview
      var joinResult = props.joinResult

      var previewNodes = []
      if (autoPreview && autoPreview.rows && autoPreview.rows.length) {
        previewNodes.push(h('div', { key: 'auto-preview', className: 'sc-samples-preview' },
          h('div', { className: 'sc-record-title' }, T('samplesAutoPreview')),
          h('div', { className: 'sc-samples-preview-list' }, autoPreview.rows.map(function (row, index) {
            return h('span', { key: 'ap-' + index, className: 'sc-samples-preview-chip' },
              row.id, ' → ', h('b', null, row.group || '—'))
          })),
          autoPreview.unmatched && autoPreview.unmatched.length
            ? h('div', { className: 'sc-hint-line' }, T('samplesAutoEmpty') + '：' + autoPreview.unmatched.slice(0, 8).join(', '))
            : null,
          h('div', { className: 'sc-inline-row' },
            h('button', { type: 'button', className: 'sc-btn sc-btn-sm sc-btn-primary', disabled: saving, onClick: props.onAutoApply }, saving ? T('samplesSaving') : T('samplesApply')),
            h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: saving, onClick: props.onAutoCancel }, T('samplesCancel')))))
      }
      if (importPreview) {
        var updateKeys = Object.keys(importPreview.updates || {}).sort()
        previewNodes.push(h('div', { key: 'import-preview', className: 'sc-samples-preview' },
          h('div', { className: 'sc-record-title' }, tf(T, 'samplesImportMapResult', formatNumber(importPreview.matched) + '/' + formatNumber(importPreview.total))),
          h('div', { className: 'sc-samples-preview-list' }, updateKeys.map(function (key) {
            return h('span', { key: 'ip-' + key, className: 'sc-samples-preview-chip' }, key, ' → ', h('b', null, importPreview.updates[key]))
          })),
          h('div', { className: 'sc-inline-row' },
            h('button', { type: 'button', className: 'sc-btn sc-btn-sm sc-btn-primary', disabled: saving, onClick: props.onImportApply }, saving ? T('samplesSaving') : T('samplesApply')),
            h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: saving, onClick: props.onImportCancel }, T('samplesCancel')))))
      }
      if (joinResult) {
        previewNodes.push(h('div', { key: 'join-result', className: 'sc-samples-preview' },
          h('div', { className: 'sc-record-title' }, tf(T, 'samplesJoinResult', (joinResult.joinedColumns || []).join(', ') || '—')),
          h('div', { className: 'sc-record-note' }, T('samplesJoinRate') + '：' + (formatMatchRate(joinResult) || '—'))))
      }

      var bodyNode = null
      if (props.loading && !rows.length) {
        bodyNode = LoadingLine({ text: T('loading') })
      } else if (props.error) {
        bodyNode = ErrorBanner({ message: props.error, onRetry: props.onReload, retryLabel: T('retry') })
      } else if (!rows.length) {
        bodyNode = h('div', { className: 'sc-record-note' }, T('samplesEmpty'))
      } else {
        bodyNode = h('div', { className: 'sc-samples-table-wrap' },
          h('table', { className: 'sc-samples-table' },
            h('thead', null, h('tr', null,
              h('th', null, T('samplesColSource')),
              h('th', { className: 'sc-td-num' }, T('samplesColCells')),
              h('th', null, T('samplesColName')),
              h('th', null, T('samplesColGroup')),
              h('th', null, T('samplesColActions')))),
            h('tbody', null, rows.map(function (row) {
              var known = false
              for (var optionIndex = 0; optionIndex < groupOptionsList.length; optionIndex++) {
                if (groupOptionsList[optionIndex].value === (row.group || '')) known = true
              }
              var optionNodes = groupOptionsList.map(function (option) {
                return h('option', { key: 'go-' + option.value, value: option.value }, option.label)
              })
              if (!known) optionNodes.push(h('option', { key: 'go-current', value: row.group || '' }, row.group || '—'))
              optionNodes.push(h('option', { key: 'go-clear', value: '__clear__' }, T('samplesClearGroup')))
              return h('tr', { key: 'sample-' + row.id },
                h('td', { className: 'sc-samples-src', title: row.id }, row.id),
                h('td', { className: 'sc-td-num' }, formatNumber(row.nCells)),
                h('td', null, h('input', {
                  className: 'sc-input sc-sample-name',
                  type: 'text',
                  value: row.name || '',
                  disabled: saving,
                  'aria-label': T('samplesColName'),
                  onChange: function (event) { props.onNameChange(row.id, event.target.value) },
                })),
                h('td', null, h('select', {
                  className: 'sc-input sc-select sc-sample-group',
                  value: row.group || '',
                  disabled: saving,
                  'aria-label': T('samplesColGroup'),
                  onChange: function (event) { props.onGroupChange(row.id, clearSampleGroupValue(event.target.value)) },
                }, optionNodes)),
                h('td', null, h('button', {
                  type: 'button',
                  className: 'sc-btn sc-btn-sm sc-sample-del',
                  disabled: saving,
                  title: T('samplesRemoveHint'),
                  onClick: function () { props.onRemoveSample(row.id) },
                }, T('samplesRemoveSample'))))
            }))))
      }

      return h('section', { className: 'sc-samples' },
        h('div', { className: 'sc-samples-head' },
          h('div', null,
            h('div', { className: 'sc-record-title' }, T('samplesTitle')),
            h('div', { className: 'sc-record-note' }, T('samplesSubtitle'))),
          h('div', { className: 'sc-samples-head-actions' },
            h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: saving || Boolean(props.loading), onClick: props.onReload },
              T('samplesReload')),
            h('button', { type: 'button', className: 'sc-btn sc-btn-sm sc-btn-primary', disabled: saving || !rows.length, onClick: props.onSave },
              saving ? T('samplesSaving') : T('samplesSave')))),
        props.notice ? h('div', { className: 'sc-samples-notice' }, h('span', null, props.notice)) : null,
        props.statusText ? h('div', { className: 'sc-samples-status' }, T('samplesStatusPrefix') + props.statusText) : null,
        (props.groups && props.groups.length)
          ? h('div', { className: 'sc-groups-row' },
              h('span', { className: 'sc-samples-tool-label' }, T('samplesGroupsTitle')),
              props.groups.map(function (group) {
                var used = isGroupUsedBySamples(props.rows, group)
                var color = (props.groupColors || {})[group]
                return h('span', { key: 'sg-' + group, className: 'sc-group-chip' },
                  color ? h('span', { className: 'sc-group-chip-color', style: { background: color } }) : null,
                  group,
                  h('button', {
                    type: 'button',
                    className: 'sc-group-chip-del',
                    disabled: used,
                    title: used ? T('samplesGroupUsed') : T('samplesDeleteGroup'),
                    onClick: function () { props.onDeleteGroup(group) },
                  }, '×'))
              }))
          : null,
        props.dirty ? h('div', { className: 'sc-hint-line' }, T('samplesDirty')) : null,
        bodyNode,
        previewNodes.length ? h('div', { className: 'sc-samples-tools' }, previewNodes) : null,
        h('div', { className: 'sc-samples-tools' },
          h('div', { className: 'sc-samples-tool' },
            h('span', { className: 'sc-samples-tool-label' }, T('samplesNewGroup')),
            h('input', {
              className: 'sc-input sc-samples-tool-input sc-samples-key',
              type: 'text',
              value: props.groupDraft || '',
              placeholder: T('samplesAddGroupPlaceholder'),
              disabled: saving,
              onChange: function (event) { props.onGroupDraft(event.target.value) },
              onKeyDown: function (event) {
                if (event.key === 'Enter') {
                  event.preventDefault()
                  props.onAddGroup()
                }
              },
            }),
            h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: saving || !String(props.groupDraft || '').trim(), onClick: props.onAddGroup }, T('samplesAdd'))),
          h('div', { className: 'sc-samples-tool' },
            h('span', { className: 'sc-samples-tool-label' }, T('samplesAutoGroup')),
            h('input', {
              className: 'sc-input sc-samples-tool-input',
              type: 'text',
              value: props.autoPattern || '',
              placeholder: T('samplesAutoPatternPlaceholder'),
              disabled: saving || Boolean(props.autoBusy),
              onChange: function (event) { props.onAutoPattern(event.target.value) },
              onKeyDown: function (event) {
                if (event.key === 'Enter') {
                  event.preventDefault()
                  props.onAutoRun()
                }
              },
            }),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: saving || Boolean(props.autoBusy) || !String(props.autoPattern || '').trim(),
              onClick: props.onAutoRun,
            }, props.autoBusy ? T('running') : T('samplesAutoRun'))),
          h('div', { className: 'sc-samples-tool' },
            h('span', { className: 'sc-samples-tool-label' }, T('samplesImportMap')),
            h('input', {
              className: 'sc-input sc-samples-tool-input',
              type: 'text',
              value: props.importPath || '',
              placeholder: T('samplesImportMapPath'),
              disabled: saving || Boolean(props.mapBusy),
              onChange: function (event) { props.onImportPath(event.target.value) },
            }),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: saving || Boolean(props.mapBusy) || !String(props.importPath || '').trim(),
              onClick: props.onImportRun,
            }, props.mapBusy ? T('running') : T('samplesImportMapRun'))),
          h('div', { className: 'sc-samples-tool' },
            h('span', { className: 'sc-samples-tool-label' }, T('samplesJoin')),
            h('input', {
              className: 'sc-input sc-samples-tool-input',
              type: 'text',
              value: props.joinPath || '',
              placeholder: T('samplesJoinPath'),
              disabled: saving || Boolean(props.joinBusy),
              onChange: function (event) { props.onJoinPath(event.target.value) },
            }),
            h('input', {
              className: 'sc-input sc-samples-tool-input sc-samples-key',
              type: 'text',
              value: props.joinKey || '',
              placeholder: T('samplesJoinKey'),
              disabled: saving || Boolean(props.joinBusy),
              onChange: function (event) { props.onJoinKey(event.target.value) },
            }),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: saving || Boolean(props.joinBusy) || !String(props.joinPath || '').trim(),
              onClick: props.onJoinRun,
            }, props.joinBusy ? T('running') : T('samplesJoinRun')))))

    }


    // ==========================================================================
    // 记录 tab（v2）：历史 / 自检 / 自动注释 / 分辨率扫描 / 导出
    // ==========================================================================

    function RecordTab(props) {
      if (!props.activeId) {
        return h('div', { className: 'sc-guide' },
          h('div', { className: 'sc-guide-hint' }, T('recordNeedDataset')))
      }

      // 1) 数据集历史
      var history = (props.history || []).slice(-MAX_HISTORY_ROWS).reverse()
      var historyNode = history.length
        ? h('div', { className: 'sc-history' }, [
            h('div', { key: 'h-head', className: 'sc-history-row sc-history-head' },
              h('span', null, T('colTime')),
              h('span', null, T('colStep')),
              h('span', null, T('colParams')),
              h('span', null, T('colSummary')),
              h('span', { className: 'sc-history-ms' }, T('colMs'))),
          ].concat(history.map(function (entry, index) {
            var view = formatHistoryEntry(entry)
            return h('div', { key: 'h-' + index, className: 'sc-history-row', title: view.title },
              h('span', { className: 'sc-history-time' }, view.time || '—'),
              h('span', { className: 'sc-history-step' }, view.title || view.step || '—'),
              h('span', { className: 'sc-history-params' }, view.paramsText || '—'),
              h('span', { className: 'sc-history-summary' }, view.summary || (view.ok ? '' : T('resultFail'))),
              h('span', { className: 'sc-history-ms' },
                (view.duration || '') + (view.ok ? '' : ' · ' + T('resultFail'))))
          })))
        : h('div', { className: 'sc-record-note' }, T('historyEmpty'))

      // 2) 自检
      var auditNode = null
      if (props.auditLoading) {
        auditNode = LoadingLine({ text: T('auditRunning') })
      } else if (props.auditError) {
        auditNode = ErrorBanner({ message: props.auditError, onRetry: props.onAudit, retryLabel: T('retry') })
      } else if (!props.audit) {
        auditNode = h('div', { className: 'sc-record-note' }, T('auditEmpty'))
      } else {
        var auditStyle = AUDIT_LEVEL_STYLE[props.audit.level] || AUDIT_LEVEL_STYLE.ok
        auditNode = h(Fragment, null,
          h('div', { className: 'sc-audit-head' },
            h('span', { className: 'sc-audit-score sc-audit-score-' + auditStyle.tone },
              tf(T, 'auditScore', props.audit.score)),
            h('span', { className: 'sc-muted' }, T(auditStyle.labelKey))),
          h('div', { className: 'sc-check-list' }, (props.audit.checks || []).map(function (check) {
            var style = AUDIT_LEVEL_STYLE[check.level] || AUDIT_LEVEL_STYLE.ok
            return h('div', { key: check.id, className: 'sc-check-row sc-check-row-' + style.tone },
              h('span', { className: 'sc-check-icon sc-check-icon-' + style.tone, 'aria-hidden': 'true' }, style.icon),
              h('div', null,
                h('div', { className: 'sc-check-title' }, check.title || check.id),
                check.detail ? h('div', { className: 'sc-check-detail' }, check.detail) : null,
                check.hint ? h('div', { className: 'sc-check-hint' }, check.hint) : null))
          })))
      }

      // 3) 自动注释建议
      var annNode = null
      if (!props.clusterDone) {
        annNode = h('div', { className: 'sc-record-note' }, T('adoptNeedCluster'))
      } else if (props.autoannLoading) {
        annNode = LoadingLine({ text: T('autoannotateRunning') })
      } else if (props.autoannError) {
        annNode = ErrorBanner({ message: props.autoannError, onRetry: props.onAutoannotate, retryLabel: T('retry') })
      } else if (!props.autoann || !(props.autoann.suggestions || []).length) {
        annNode = h('div', { className: 'sc-record-note' }, T('autoannotateEmpty'))
      } else {
        var suggestions = props.autoann.suggestions || []
        annNode = h(Fragment, null,
          h('div', { className: 'sc-toolbar' },
            h('span', { className: 'sc-muted' }, T('autoannotatePanel') + ' · ' + (props.autoann.panel || '—')),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: props.autoannLoading,
              onClick: props.onAutoannotate,
            }, T('runAutoannotate')),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm sc-btn-primary',
              disabled: props.adopting,
              onClick: props.onAdoptAll,
            }, props.adopting ? T('running') : T('adoptAll'))),
          h('div', { className: 'sc-ann-list' }, suggestions.map(function (suggestion, index) {
            var adopted = props.adopted && props.adopted[suggestion.cluster]
            var evidence = (suggestion.evidence || []).slice(0, 5)
            return h('div', { key: 'ann-' + suggestion.cluster + '-' + index, className: 'sc-ann-row' },
              h('span', { className: 'sc-ann-cluster' }, suggestion.cluster),
              h('span', { className: 'sc-ann-label', title: suggestion.label }, suggestion.label || '—'),
              h('div', { className: 'sc-ann-meta' },
                h('div', { className: 'sc-ann-evidence' },
                  T('autoannotateScore') + ' ' + formatNumber(suggestion.score, 2) +
                  ' · ' + T('autoannotateMargin') + ' ' + formatNumber(suggestion.margin, 2) +
                  (suggestion.alternative ? ' · ' + T('autoannotateAlternative') + ' ' + suggestion.alternative : '')),
                evidence.length
                  ? h('div', { className: 'sc-ann-evidence' }, T('autoannotateEvidence') + ' · ' + evidence.join(', '))
                  : null),
              h('div', { className: 'sc-ann-actions' },
                h('button', {
                  type: 'button',
                  className: 'sc-btn sc-btn-sm',
                  disabled: adopted || props.adopting,
                  onClick: function () { props.onAdoptCluster(suggestion) },
                }, adopted ? T('adopted') : T('adopt'))))
          })))
        annNode = h(Fragment, null, annNode,
          props.adoptError ? h('div', { className: 'sc-inline-err' }, props.adoptError) : null)
      }

      // 4) 分辨率扫描
      var scanNode = null
      if (props.scanLoading) {
        scanNode = LoadingLine({ text: T('scanRunning') })
      } else if (props.scanError) {
        scanNode = ErrorBanner({ message: props.scanError, onRetry: props.onScan, retryLabel: T('retry') })
      } else if (!props.scan) {
        scanNode = h('div', { className: 'sc-record-note' }, T('scanEmpty'))
      } else {
        var currentText = props.scan.current && props.scan.current.resolution !== null
          ? T('scanCurrent') + '：resolution=' + formatNumber(props.scan.current.resolution, 2) + ' · ' +
            T('colClusters') + ' ' + formatNumber(props.scan.current.nClusters)
          : ''
        scanNode = h(Fragment, null,
          currentText ? h('div', { className: 'sc-muted' }, currentText) : null,
          h('div', { className: 'sc-table-wrap' },
            h('table', { className: 'sc-table' },
              h('thead', null, h('tr', null,
                h('th', null, T('colResolution')),
                h('th', { className: 'sc-td-num' }, T('colClusters')),
                h('th', { className: 'sc-td-num' }, T('colAri')))),
              h('tbody', null, (props.scan.results || []).map(function (row, index) {
                return h('tr', { key: 'scan-' + index },
                  h('td', null, formatNumber(row.resolution, 2)),
                  h('td', { className: 'sc-td-num' }, formatNumber(row.nClusters)),
                  h('td', { className: 'sc-td-num' }, row.ari === null || row.ari === undefined ? '—' : row.ari.toFixed(3)))
              })))))
      }

      // 5) 导出
      var exportNode = h(Fragment, null,
        h('div', { className: 'sc-export-btns' }, EXPORT_KINDS.map(function (kind) {
          var busy = props.exportLoadingKind === kind.id
          return h('button', {
            key: kind.id,
            type: 'button',
            className: 'sc-btn sc-btn-sm',
            title: T(kind.labelKey) + ' · ' + kind.hint,
            disabled: Boolean(props.exportLoadingKind),
            onClick: function () { props.onExport(kind.id) },
          }, busy ? T('exportRunning') : T(kind.labelKey))
        })),
        props.exportError ? h('div', { className: 'sc-inline-err' }, props.exportError) : null,
        props.exportResult && props.exportResult.file
          ? h('div', { className: 'sc-export-result' },
              h('div', { className: 'sc-muted' }, T('exportDone') + (props.exportResult.what || '')),
              h('div', { className: 'sc-path' }, props.exportResult.file),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm',
                onClick: props.onCopyPath,
              }, props.copiedKey === 'export-path' ? T('copied') : T('copyPath')))
          : null)

      return h('div', { className: 'sc-record' },
        props.samplesPanel
          ? h('div', { className: 'sc-record-block' }, props.samplesPanel)
          : null,
        h('div', { className: 'sc-record-block' },
          h('div', { className: 'sc-record-head' },
            h('span', { className: 'sc-record-title' }, T('historyTitle'))),
          historyNode),
        h('div', { className: 'sc-record-block' },
          h('div', { className: 'sc-record-head' },
            h('span', { className: 'sc-record-title' }, T('auditTitle')),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: props.auditLoading,
              onClick: props.onAudit,
            }, props.auditLoading ? T('auditRunning') : T('runAudit'))),
          auditNode),
        h('div', { className: 'sc-record-block' },
          h('div', { className: 'sc-record-head' },
            h('span', { className: 'sc-record-title' }, T('autoannotateTitle')),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: props.autoannLoading || !props.clusterDone,
              onClick: props.onAutoannotate,
            }, props.autoannLoading ? T('autoannotateRunning') : T('runAutoannotate'))),
          annNode),
        h('div', { className: 'sc-record-block' },
          h('div', { className: 'sc-record-head' },
            h('span', { className: 'sc-record-title' }, T('resolutionTitle'))),
          h('div', { className: 'sc-scan-form' },
            h('input', {
              className: 'sc-input sc-scan-input',
              type: 'text',
              value: props.scanInput,
              placeholder: T('resolutionPlaceholder'),
              onChange: function (event) { props.onScanInputChange(event.target.value) },
              onKeyDown: function (event) {
                if (event.key === 'Enter') {
                  event.preventDefault()
                  props.onScan()
                }
              },
            }),
            h('button', {
              type: 'button',
              className: 'sc-btn sc-btn-sm',
              disabled: props.scanLoading,
              onClick: props.onScan,
            }, props.scanLoading ? T('scanRunning') : T('runScan'))),
          h('div', { className: 'sc-record-note' }, T('scanHint')),
          scanNode),
        h('div', { className: 'sc-record-block' },
          h('div', { className: 'sc-record-head' },
            h('span', { className: 'sc-record-title' }, T('exportTitle'))),
          exportNode))
    }

    // ==========================================================================
    // 日志 tab
    // ==========================================================================

    function LogTab(props) {
      var rows = (props.entries || EMPTY_LIST).map(function (entry, index) {
        return h('div', { key: 'log-' + index, className: 'sc-log-row' },
          h('span', { className: 'sc-log-time' }, formatClock(entry.at) || '—'),
          h('span', { className: 'sc-log-actor', 'aria-hidden': 'true' }, activityIcon(entry.actor)),
          h('span', { className: 'sc-log-action', title: entry.action }, entry.action || '—'),
          h('span', { className: 'sc-log-detail', title: entry.detail }, entry.detail || ''),
          h('span', { className: 'sc-log-ms' },
            (entry.ok ? '' : '✕ ') + (entry.ms !== null && entry.ms !== undefined ? formatDuration(entry.ms) : '')))
      })
      return h(Fragment, null,
        h('div', { className: 'sc-toolbar' },
          h('span', { className: 'sc-muted' }, T('activity')),
          h('button', { type: 'button', className: 'sc-btn sc-btn-sm', disabled: props.loading, onClick: props.onRefresh },
            T('refresh'))),
        props.error
          ? ErrorBanner({ message: props.error, onRetry: props.onRefresh, retryLabel: T('retry') })
          : null,
        rows.length ? h('div', { className: 'sc-log-list' }, rows) : Notice({ text: T('logEmpty') }))
    }

    // ==========================================================================
    // 工作台主体：状态机 + 轮询 + 所有 API 动作
    // ==========================================================================

    var clientApi = createApi()

    function WorkbenchBody() {
      // ---- 基础数据 ----
      var stateRes = React.useState(null)
      var state = stateRes[0]
      var setState = stateRes[1]
      var healthEngineRes = React.useState(null)
      var healthEngine = healthEngineRes[0]
      var setHealthEngine = healthEngineRes[1]
      var bootingRes = React.useState(true)
      var booting = bootingRes[0]
      var setBooting = bootingRes[1]
      var bootErrorRes = React.useState('')
      var bootError = bootErrorRes[0]
      var setBootError = bootErrorRes[1]
      var pollErrorRes = React.useState('')
      var pollError = pollErrorRes[0]
      var setPollError = pollErrorRes[1]
      var tabRes = React.useState('embedding')
      var tab = tabRes[0]
      var setTab = tabRes[1]
      var busyRes = React.useState({})
      var busy = busyRes[0]
      var setBusy = busyRes[1]
      var errorsRes = React.useState({})
      var errors = errorsRes[0]
      var setErrors = errorsRes[1]
      var previewRes = React.useState(null)
      var preview = previewRes[0]
      var setPreview = previewRes[1]

      // ---- 流程卡片 ----
      var summariesRes = React.useState({})
      var stepSummaries = summariesRes[0]
      var setStepSummaries = summariesRes[1]

      // ---- 嵌入图 ----
      var colorByRes = React.useState({ kind: 'cluster', key: 'clusters' })
      var colorBy = colorByRes[0]
      var setColorBy = colorByRes[1]
      var geneQueryRes = React.useState('')
      var geneQuery = geneQueryRes[0]
      var setGeneQuery = geneQueryRes[1]
      var embeddingRes = React.useState(null)
      var embedding = embeddingRes[0]
      var setEmbedding = embeddingRes[1]
      var embeddingLoadingRes = React.useState(false)
      var embeddingLoading = embeddingLoadingRes[0]
      var setEmbeddingLoading = embeddingLoadingRes[1]
      var embeddingErrorRes = React.useState('')
      var embeddingError = embeddingErrorRes[0]
      var setEmbeddingError = embeddingErrorRes[1]
      var highlightRes = React.useState(null)
      var highlight = highlightRes[0]
      var setHighlight = highlightRes[1]
      var embedNonceRes = React.useState(0)
      var embedNonce = embedNonceRes[0]
      var setEmbedNonce = embedNonceRes[1]

      // ---- QC ----
      var qcMetricsRes = React.useState(['n_genes', 'n_counts', 'pct_counts_mt'])
      var qcMetrics = qcMetricsRes[0]
      var setQcMetrics = qcMetricsRes[1]
      var qcFigureRes = React.useState(null)
      var qcFigure = qcFigureRes[0]
      var setQcFigure = qcFigureRes[1]

      // ---- 标记基因 / 聚类 ----
      var markersRes = React.useState(null)
      var markers = markersRes[0]
      var setMarkers = markersRes[1]
      var markersLoadingRes = React.useState(false)
      var markersLoading = markersLoadingRes[0]
      var setMarkersLoading = markersLoadingRes[1]
      var markersErrorRes = React.useState('')
      var markersError = markersErrorRes[0]
      var setMarkersError = markersErrorRes[1]
      var markersTopRes = React.useState(10)
      var markersTop = markersTopRes[0]
      var setMarkersTop = markersTopRes[1]
      var markersMethodRes = React.useState('t-test')
      var markersMethod = markersMethodRes[0]
      var setMarkersMethod = markersMethodRes[1]
      var markersNonceRes = React.useState(0)
      var markersNonce = markersNonceRes[0]
      var setMarkersNonce = markersNonceRes[1]
      var dotplotRes = React.useState(null)
      var dotplot = dotplotRes[0]
      var setDotplot = dotplotRes[1]
      var clusterBarRes = React.useState(null)
      var clusterBar = clusterBarRes[0]
      var setClusterBar = clusterBarRes[1]
      var clusterBarLoadingRes = React.useState(false)
      var clusterBarLoading = clusterBarLoadingRes[0]
      var setClusterBarLoading = clusterBarLoadingRes[1]
      var clusterBarErrorRes = React.useState('')
      var clusterBarError = clusterBarErrorRes[0]
      var setClusterBarError = clusterBarErrorRes[1]
      var clusterBarNonceRes = React.useState(0)
      var clusterBarNonce = clusterBarNonceRes[0]
      var setClusterBarNonce = clusterBarNonceRes[1]
      var annotationRes = React.useState({})
      var annotation = annotationRes[0]
      var setAnnotation = annotationRes[1]
      var annotationSavedRes = React.useState('')
      var annotationSaved = annotationSavedRes[0]
      var setAnnotationSaved = annotationSavedRes[1]
      var adoptedRes = React.useState({})
      var adopted = adoptedRes[0]
      var setAdopted = adoptedRes[1]

      // ---- 日志 ----
      var logRes = React.useState(null)
      var log = logRes[0]
      var setLog = logRes[1]
      var logLoadingRes = React.useState(false)
      var logLoading = logLoadingRes[0]
      var setLogLoading = logLoadingRes[1]
      var logErrorRes = React.useState('')
      var logError = logErrorRes[0]
      var setLogError = logErrorRes[1]
      var logNonceRes = React.useState(0)
      var logNonce = logNonceRes[0]
      var setLogNonce = logNonceRes[1]

      // ---- 记录 tab（v2）----
      var auditRes = React.useState(null)
      var audit = auditRes[0]
      var setAudit = auditRes[1]
      var auditErrorRes = React.useState('')
      var auditError = auditErrorRes[0]
      var setAuditError = auditErrorRes[1]
      var autoannRes = React.useState(null)
      var autoann = autoannRes[0]
      var setAutoann = autoannRes[1]
      var autoannErrorRes = React.useState('')
      var autoannError = autoannErrorRes[0]
      var setAutoannError = autoannErrorRes[1]
      var scanRes = React.useState(null)
      var scan = scanRes[0]
      var setScan = scanRes[1]
      var scanErrorRes = React.useState('')
      var scanError = scanErrorRes[0]
      var setScanError = scanErrorRes[1]
      var scanInputRes = React.useState('0.4,0.8,1.2')
      var scanInput = scanInputRes[0]
      var setScanInput = scanInputRes[1]
      var exportResultRes = React.useState(null)
      var exportResult = exportResultRes[0]
      var setExportResult = exportResultRes[1]
      var exportErrorRes = React.useState('')
      var exportError = exportErrorRes[0]
      var setExportError = exportErrorRes[1]
      var exportWhatRes = React.useState('')
      var exportWhat = exportWhatRes[0]
      var setExportWhat = exportWhatRes[1]

      // ---- 预设流程 ----
      var recipeIdRes = React.useState('')
      var recipeId = recipeIdRes[0]
      var setRecipeId = recipeIdRes[1]
      var recipeResultRes = React.useState(null)
      var recipeResult = recipeResultRes[0]
      var setRecipeResult = recipeResultRes[1]

      // ---- 复制提示 ----
      var copiedKeyRes = React.useState('')
      var copiedKey = copiedKeyRes[0]
      var setCopiedKey = copiedKeyRes[1]
      var copyErrorRes = React.useState('')
      var copyError = copyErrorRes[0]
      var setCopyError = copyErrorRes[1]
      var copyTimerRef = React.useRef(null)
      var healthRequestedRef = React.useRef(false)

      // ---- v3：导入数据面板 ----
      var dialogOpenRes = React.useState(false)
      var dialogOpen = dialogOpenRes[0]
      var setDialogOpen = dialogOpenRes[1]
      var dialogNonceRes = React.useState(0)
      var dialogNonce = dialogNonceRes[0]
      var setDialogNonce = dialogNonceRes[1]
      var importBusyRes = React.useState(false)
      var importBusy = importBusyRes[0]
      var setImportBusy = importBusyRes[1]
      var importErrorRes = React.useState('')
      var importError = importErrorRes[0]
      var setImportError = importErrorRes[1]
      var importCandidatesRes = React.useState(EMPTY_LIST)
      var importCandidates = importCandidatesRes[0]
      var setImportCandidates = importCandidatesRes[1]
      var importCandidatesDisabledRes = React.useState(false)
      var importCandidatesDisabled = importCandidatesDisabledRes[0]
      var setImportCandidatesDisabled = importCandidatesDisabledRes[1]
      var lastImportRef = React.useRef(null)
      var formatsRes = React.useState(null)
      var formats = formatsRes[0]
      var setFormats = formatsRes[1]
      var formatsLoadingRes = React.useState(false)
      var formatsLoading = formatsLoadingRes[0]
      var setFormatsLoading = formatsLoadingRes[1]
      var formatsErrorRes = React.useState('')
      var formatsError = formatsErrorRes[0]
      var setFormatsError = formatsErrorRes[1]
      var inspectBusyRes = React.useState(false)
      var inspectBusy = inspectBusyRes[0]
      var setInspectBusy = inspectBusyRes[1]
      var inspectResultRes = React.useState(null)
      var inspectResult = inspectResultRes[0]
      var setInspectResult = inspectResultRes[1]
      var inspectErrorRes = React.useState('')
      var inspectError = inspectErrorRes[0]
      var setInspectError = inspectErrorRes[1]

      // ---- v3：高级分析 ----
      var analysesPayloadRes = React.useState(null)
      var analysesPayload = analysesPayloadRes[0]
      var setAnalysesPayload = analysesPayloadRes[1]
      var analysesLoadingRes = React.useState(false)
      var analysesLoading = analysesLoadingRes[0]
      var setAnalysesLoading = analysesLoadingRes[1]
      var analysesErrorRes = React.useState('')
      var analysesError = analysesErrorRes[0]
      var setAnalysesError = analysesErrorRes[1]
      var analysesNonceRes = React.useState(0)
      var analysesNonce = analysesNonceRes[0]
      var setAnalysesNonce = analysesNonceRes[1]
      var analysesKeyRef = React.useRef('')
      var selectedAnalysisRes = React.useState(null)
      var selectedAnalysis = selectedAnalysisRes[0]
      var setSelectedAnalysis = selectedAnalysisRes[1]
      var analysisValuesRes = React.useState({})
      var analysisValues = analysisValuesRes[0]
      var setAnalysisValues = analysisValuesRes[1]
      var analysisResultRes = React.useState(null)
      var analysisResult = analysisResultRes[0]
      var setAnalysisResult = analysisResultRes[1]
      var analysisErrorRes = React.useState('')
      var analysisError = analysisErrorRes[0]
      var setAnalysisError = analysisErrorRes[1]
      var previewNonceRes = React.useState(0)
      var previewNonce = previewNonceRes[0]
      var setPreviewNonce = previewNonceRes[1]

      // ---- v4：样本 & 分组 ----
      var obsPayloadRes = React.useState(null)
      var obsPayload = obsPayloadRes[0]
      var setObsPayload = obsPayloadRes[1]
      var obsLoadingRes = React.useState(false)
      var obsLoading = obsLoadingRes[0]
      var setObsLoading = obsLoadingRes[1]
      var obsErrorRes = React.useState('')
      var obsError = obsErrorRes[0]
      var setObsError = obsErrorRes[1]
      var samplesNonceRes = React.useState(0)
      var samplesNonce = samplesNonceRes[0]
      var setSamplesNonce = samplesNonceRes[1]
      var sampleRowsStateRes = React.useState(EMPTY_LIST)
      var sampleRowsState = sampleRowsStateRes[0]
      var setSampleRows = sampleRowsStateRes[1]
      var sampleGroupsRes = React.useState(EMPTY_LIST)
      var sampleGroups = sampleGroupsRes[0]
      var setSampleGroups = sampleGroupsRes[1]
      var sampleDirtyRes = React.useState(false)
      var sampleDirty = sampleDirtyRes[0]
      var setSampleDirty = sampleDirtyRes[1]
      var samplesMessageRes = React.useState('')
      var samplesMessage = samplesMessageRes[0]
      var setSamplesMessage = samplesMessageRes[1]
      var groupDraftRes = React.useState('')
      var groupDraft = groupDraftRes[0]
      var setGroupDraft = groupDraftRes[1]
      var autoPatternRes = React.useState('')
      var autoPattern = autoPatternRes[0]
      var setAutoPattern = autoPatternRes[1]
      var autoPreviewRes = React.useState(null)
      var autoPreview = autoPreviewRes[0]
      var setAutoPreview = autoPreviewRes[1]
      var sampleMapPathRes = React.useState('')
      var sampleMapPath = sampleMapPathRes[0]
      var setSampleMapPath = sampleMapPathRes[1]
      var sampleMapPreviewRes = React.useState(null)
      var sampleMapPreview = sampleMapPreviewRes[0]
      var setSampleMapPreview = sampleMapPreviewRes[1]
      var joinPathRes = React.useState('')
      var joinPath = joinPathRes[0]
      var setJoinPath = joinPathRes[1]
      var joinKeyRes = React.useState('')
      var joinKey = joinKeyRes[0]
      var setJoinKey = joinKeyRes[1]
      var joinResultRes = React.useState(null)
      var joinResult = joinResultRes[0]
      var setJoinResult = joinResultRes[1]
      var mergeNoticeRes = React.useState(null)
      var mergeNotice = mergeNoticeRes[0]
      var setMergeNotice = mergeNoticeRes[1]

      // ---- v5：数据根目录扫描 → 样本表 ----
      var scanRootRes = React.useState('')
      var scanRoot = scanRootRes[0]
      var setScanRoot = scanRootRes[1]
      var scanRecursiveRes = React.useState(true)
      var scanRecursive = scanRecursiveRes[0]
      var setScanRecursive = scanRecursiveRes[1]
      var scanLoadingRes = React.useState(false)
      var scanLoading = scanLoadingRes[0]
      var setScanLoading = scanLoadingRes[1]
      var scanErrorRes = React.useState('')
      var scanError = scanErrorRes[0]
      var setScanError = scanErrorRes[1]
      var scanResultRes = React.useState(null)
      var scanResult = scanResultRes[0]
      var setScanResult = scanResultRes[1]
      var scanRowsRes = React.useState(EMPTY_LIST)
      var scanRows = scanRowsRes[0]
      var setScanRows = scanRowsRes[1]
      var scanGroupsRes = React.useState(EMPTY_LIST)
      var scanGroups = scanGroupsRes[0]
      var setScanGroups = scanGroupsRes[1]
      var scanGroupColorsRes = React.useState({})
      var scanGroupColors = scanGroupColorsRes[0]
      var setScanGroupColors = scanGroupColorsRes[1]
      var scanNewGroupRes = React.useState('')
      var scanNewGroup = scanNewGroupRes[0]
      var setScanNewGroup = scanNewGroupRes[1]
      var scanNewGroupColorRes = React.useState('#4e79a7')
      var scanNewGroupColor = scanNewGroupColorRes[0]
      var setScanNewGroupColor = scanNewGroupColorRes[1]
      var scanBatchGroupRes = React.useState('')
      var scanBatchGroup = scanBatchGroupRes[0]
      var setScanBatchGroup = scanBatchGroupRes[1]
      var scanConfirmErrorRes = React.useState('')
      var scanConfirmError = scanConfirmErrorRes[0]
      var setScanConfirmError = scanConfirmErrorRes[1]
      var geoGseRes = React.useState('')
      var geoGse = geoGseRes[0]
      var setGeoGse = geoGseRes[1]
      var geoLoadingRes = React.useState(false)
      var geoLoading = geoLoadingRes[0]
      var setGeoLoading = geoLoadingRes[1]
      var geoResultRes = React.useState(null)
      var geoResult = geoResultRes[0]
      var setGeoResult = geoResultRes[1]
      var geoErrorRes = React.useState('')
      var geoError = geoErrorRes[0]
      var setGeoError = geoErrorRes[1]
      var openedImportStateRes = React.useState(null)
      var openedImportState = openedImportStateRes[0]
      var setOpenedImportState = openedImportStateRes[1]
      var importStateDismissedRes = React.useState({})
      var importStateDismissed = importStateDismissedRes[0]
      var setImportStateDismissed = importStateDismissedRes[1]

      // ---- v6：降维聚类 / 手动注释 / Agent 输入 ----
      var dimEngineRes = React.useState('seurat')
      var dimEngine = dimEngineRes[0]
      var setDimEngine = dimEngineRes[1]
      var dimHarmonyRes = React.useState(true)
      var dimHarmony = dimHarmonyRes[0]
      var setDimHarmony = dimHarmonyRes[1]
      var dimBatchKeyRes = React.useState('')
      var dimBatchKey = dimBatchKeyRes[0]
      var setDimBatchKey = dimBatchKeyRes[1]
      var dimResolutionRes = React.useState(DIMRED_DEFAULTS.resolution)
      var dimResolution = dimResolutionRes[0]
      var setDimResolution = dimResolutionRes[1]
      var dimNpcsRes = React.useState(DIMRED_DEFAULTS.nPcs)
      var dimNpcs = dimNpcsRes[0]
      var setDimNpcs = dimNpcsRes[1]
      var dimNeighborsRes = React.useState(DIMRED_DEFAULTS.neighbors)
      var dimNeighbors = dimNeighborsRes[0]
      var setDimNeighbors = dimNeighborsRes[1]
      var dimTopGenesRes = React.useState(DIMRED_DEFAULTS.nTopGenes)
      var dimTopGenes = dimTopGenesRes[0]
      var setDimTopGenes = dimTopGenesRes[1]
      var dimWarningsRes = React.useState(EMPTY_LIST)
      var dimWarnings = dimWarningsRes[0]
      var setDimWarnings = dimWarningsRes[1]
      var selectedClusterRes = React.useState(null)
      var selectedCluster = selectedClusterRes[0]
      var setSelectedCluster = selectedClusterRes[1]
      var selectedGeneRes = React.useState('')
      var selectedGene = selectedGeneRes[0]
      var setSelectedGene = selectedGeneRes[1]
      var annotPlotFigureRes = React.useState(null)
      var annotPlotFigure = annotPlotFigureRes[0]
      var setAnnotPlotFigure = annotPlotFigureRes[1]
      var annotCsvOpenRes = React.useState(false)
      var annotCsvOpen = annotCsvOpenRes[0]
      var setAnnotCsvOpen = annotCsvOpenRes[1]
      var annotCsvTextRes = React.useState('')
      var annotCsvText = annotCsvTextRes[0]
      var setAnnotCsvText = annotCsvTextRes[1]
      var agentInputRes = React.useState('')
      var agentInput = agentInputRes[0]
      var setAgentInput = agentInputRes[1]
      var agentSuggestionsRes = React.useState(EMPTY_LIST)
      var agentSuggestions = agentSuggestionsRes[0]
      var setAgentSuggestions = agentSuggestionsRes[1]
      var agentResultsRes = React.useState({})
      var agentResults = agentResultsRes[0]
      var setAgentResults = agentResultsRes[1]
      var agentBusyIdRes = React.useState('')
      var agentBusyId = agentBusyIdRes[0]
      var setAgentBusyId = agentBusyIdRes[1]
      var agentFallbackRes = React.useState(false)
      var agentFallback = agentFallbackRes[0]
      var setAgentFallback = agentFallbackRes[1]

      var normalized = state || normalizeState(null)
      var activeId = normalized.active
      var activeDataset = findDataset(normalized, activeId)
      var doneSteps = activeDataset ? activeDataset.steps : EMPTY_LIST
      var stepsKey = (activeId || '') + '|' + doneSteps.join(',')
      var engine = chooseEngine(normalized.engine, healthEngine)
      var recipes = normalized.recipes || EMPTY_LIST
      var recipeProgress = normalized.recipeProgress
      var suggestions = normalized.suggestions || EMPTY_LIST
      var activity = normalized.activity || EMPTY_LIST
      var history = activeDataset ? activeDataset.history : EMPTY_LIST
      var annotationMapping = {}
      Object.keys(annotation).forEach(function (clusterId) {
        var name = String(annotation[clusterId] || '').trim()
        if (name) annotationMapping[clusterId] = name
      })
      var annotationCount = Object.keys(annotationMapping).length
      var markerGroups = markers ? markers.groups : EMPTY_LIST
      var runningStepIds = []
      Object.keys(busy).forEach(function (key) {
        if (busy[key] && key.indexOf('step:') === 0) runningStepIds.push(key.slice(5))
      })
      // 完成后的一行摘要：本次会话的 run 结果优先；刷新页面后从 dataset.history 回填
      var historySummaries = {}
      if (activeDataset && activeDataset.history.length) {
        activeDataset.history.forEach(function (entry) {
          if (!entry || entry.ok === false) return
          var historyStep = String(entry.step || '')
          if (!STEP_BY_ID[historyStep]) return
          var historyText = typeof entry.summary === 'string' ? entry.summary : ''
          if (!historyText) return
          historySummaries[historyStep] = { title: stepTitleOf(historyStep), bullets: [historyText] }
        })
      }
      var mergedSummaries = Object.assign({}, historySummaries, stepSummaries)
      var obsColumns = mergeObsColumns(
        preview ? preview.obsColumns : EMPTY_LIST,
        activeDataset ? activeDataset.obsColumns : EMPTY_LIST,
        analysesPayload ? analysesPayload.columns : EMPTY_LIST,
        obsPayload ? obsPayload.columns : EMPTY_LIST,
      )
      var obsColumnsKey = obsColumns.map(function (column) { return column.key }).join(',')
      var colorOptions = buildColorOptions(obsColumns)
      var analysesGroups = groupAnalyses(analysesPayload ? analysesPayload.analyses : EMPTY_LIST)
      var groupColumn = (obsPayload && obsPayload.groupColumn)
        || (activeDataset && activeDataset.groupColumn)
        || 'group'
      var sampleRows = sampleRowsState
      var sampleGroupOptions = groupOptions({ samples: sampleRows, groups: sampleGroups })
      var sampleStatusText = describeGroupState(sampleGroups, groupColumn)
      var importStateSource = null
      if (preview && preview.importState && preview.importState.steps.length) importStateSource = preview.importState
      else if (openedImportState && openedImportState.state && (!openedImportState.datasetId || openedImportState.datasetId === activeId)) importStateSource = openedImportState.state
      else if (activeDataset && activeDataset.importState && activeDataset.importState.steps.length) importStateSource = activeDataset.importState
      var importStateNotes = (preview && preview.notes && preview.notes.length)
        ? preview.notes
        : (activeDataset ? activeDataset.notes : EMPTY_LIST)
      var importStateHintText = importStateSource ? formatImportStateHint(importStateSource, importStateNotes) : ''
      var importStateVisible = Boolean(importStateHintText) && importStateDismissed[activeId || ''] !== true
      var dimredSettings = normalizeDimredSettings({
        engine: dimEngine,
        harmony: dimHarmony,
        batchKey: dimBatchKey,
        resolution: dimResolution,
        nPcs: dimNpcs,
        neighbors: dimNeighbors,
        nTopGenes: dimTopGenes,
      })
      var batchColumnOptions = obsColumns.map(function (column) { return column.key })

      // ---- effects ----

      // 首次加载：state + health（health 会把 Python worker 预热起来）
      React.useEffect(function () {
        var cancelled = false
        setBooting(true)
        clientApi.getState().then(function (payload) {
          if (cancelled) return
          var next = normalizeState(payload)
          setState(next)
          setBooting(false)
          if (!next.engine.ready) loadHealth()
        }).catch(function (error) {
          if (cancelled) return
          setBootError(formatError(error) || T('requestFailed'))
          setBooting(false)
          loadHealth()
        })
        return function () {
          cancelled = true
          if (copyTimerRef.current) {
            clearTimeout(copyTimerRef.current)
            copyTimerRef.current = null
          }
        }
      }, [])

      // 轮询：页面可见时每 2.5s 拉一次 state；不可见时跳过（visibilitychange 回来立刻补一次）
      React.useEffect(function () {
        var cancelled = false
        function isHidden() {
          return typeof document !== 'undefined' && document.visibilityState === 'hidden'
        }
        function tick() {
          if (cancelled || isHidden()) return
          clientApi.getState().then(function (payload) {
            if (cancelled) return
            setState(normalizeState(payload))
            setPollError('')
          }).catch(function (error) {
            if (cancelled || isAbortError(error)) return
            setPollError(formatError(error) || T('pollFailed'))
          })
        }
        var timer = setInterval(tick, POLL_MS)
        function onVisibilityChange() {
          if (!isHidden()) tick()
        }
        if (typeof document !== 'undefined' && document.addEventListener) {
          document.addEventListener('visibilitychange', onVisibilityChange)
        }
        return function () {
          cancelled = true
          clearInterval(timer)
          if (typeof document !== 'undefined' && document.removeEventListener) {
            document.removeEventListener('visibilitychange', onVisibilityChange)
          }
        }
      }, [])

      // 切换数据集：清掉所有与数据集绑定的缓存
      React.useEffect(function () {
        setPreview(null)
        setStepSummaries({})
        setErrors({})
        setEmbedding(null)
        setEmbeddingError('')
        setHighlight(null)
        setQcFigure(null)
        setMarkers(null)
        setMarkersError('')
        setDotplot(null)
        setClusterBar(null)
        setClusterBarError('')
        setAnnotation({})
        setAnnotationSaved('')
        setAdopted({})
        setLog(null)
        setLogError('')
        setAudit(null)
        setAuditError('')
        setAutoann(null)
        setAutoannError('')
        setScan(null)
        setScanError('')
        setExportResult(null)
        setExportError('')
        setCopyError('')
        setRecipeResult(null)
        setAnalysesPayload(null)
        analysesKeyRef.current = ''
        setAnalysesLoading(false)
        setAnalysesError('')
        setSelectedAnalysis(null)
        setAnalysisValues({})
        setAnalysisResult(null)
        setAnalysisError('')
        setImportError('')
        setImportCandidates(EMPTY_LIST)
        setImportCandidatesDisabled(false)
        setInspectResult(null)
        setInspectError('')
        setObsPayload(null)
        setObsLoading(false)
        setObsError('')
        setSampleRows(EMPTY_LIST)
        setSampleGroups(EMPTY_LIST)
        setSampleDirty(false)
        setSamplesMessage('')
        setGroupDraft('')
        setAutoPattern('')
        setAutoPreview(null)
        setSampleMapPath('')
        setSampleMapPreview(null)
        setJoinPath('')
        setJoinKey('')
        setJoinResult(null)
        setMergeNotice(null)
        setScanError('')
        setScanConfirmError('')
        setGeoError('')
        setScanResult(null)
        setScanRows(EMPTY_LIST)
        setScanGroups(EMPTY_LIST)
        setScanGroupColors({})
        setScanNewGroup('')
        setGeoGse('')
        setGeoResult(null)
        setOpenedImportState(null)
        setDimWarnings(EMPTY_LIST)
        setSelectedCluster(null)
        setSelectedGene('')
        setAnnotPlotFigure(null)
        setAnnotCsvOpen(false)
        setAnnotCsvText('')
        setAgentInput('')
        setAgentSuggestions(EMPTY_LIST)
        setAgentResults({})
        setAgentBusyId('')
        setAgentFallback(false)
      }, [activeId])

      // preview（基因清单 / obs 列）：steps 变了也刷新
      React.useEffect(function () {
        if (!activeId) return undefined
        var cancelled = false
        clientApi.preview(activeId).then(function (payload) {
          if (!cancelled) setPreview(normalizePreview(payload))
        }).catch(function () {
          if (!cancelled) setPreview(null)
        })
        return function () { cancelled = true }
      }, [activeId, stepsKey, previewNonce])

      // v4：记录 tab 打开时拉样本/分组表（GET /obs）
      React.useEffect(function () {
        if (tab !== 'record' || !activeId) return undefined
        var cancelled = false
        setObsLoading(true)
        setObsError('')
        clientApi.getObs(activeId).then(function (payload) {
          if (cancelled) return
          var normalizedObs = normalizeObsPayload(payload)
          setObsPayload(normalizedObs)
          setSampleRows(parseSampleRows(normalizedObs))
          setSampleGroups(normalizedObs.groups.slice())
          setSampleDirty(false)
          setObsLoading(false)
        }).catch(function (error) {
          if (cancelled || isAbortError(error)) return
          setObsError(formatError(error) || T('requestFailed'))
          setObsLoading(false)
        })
        return function () { cancelled = true }
      }, [tab, activeId, stepsKey, samplesNonce])

      // v3：导入面板的格式清单（打开面板时按需拉；失败回落内置清单）
      React.useEffect(function () {
        if (!dialogOpen || formats) return undefined
        var cancelled = false
        setFormatsLoading(true)
        setFormatsError('')
        clientApi.getFormats().then(function (payload) {
          if (cancelled) return
          setFormats(normalizeFormats(payload))
          setFormatsLoading(false)
        }).catch(function (error) {
          if (cancelled || isAbortError(error)) return
          setFormats(normalizeFormats(null))
          setFormatsError(formatError(error) || T('requestFailed'))
          setFormatsLoading(false)
        })
        return function () { cancelled = true }
      }, [dialogOpen, formats])

      // v3：高级分析目录（tab 打开 / 数据集或 obs 列变化 / 手动刷新时拉）
      React.useEffect(function () {
        if (tab !== 'analysis' || !activeId) return undefined
        var key = activeId + '|' + stepsKey + '|' + obsColumnsKey + '|' + analysesNonce
        if (analysesPayload && analysesKeyRef.current === key) return undefined
        var cancelled = false
        setAnalysesLoading(true)
        setAnalysesError('')
        clientApi.getAnalyses(activeId).then(function (payload) {
          if (cancelled) return
          analysesKeyRef.current = key
          setAnalysesPayload(normalizeAnalysesPayload(payload))
          setAnalysesLoading(false)
        }).catch(function (error) {
          if (cancelled || isAbortError(error)) return
          setAnalysesError(formatError(error) || T('requestFailed'))
          setAnalysesLoading(false)
        })
        return function () { cancelled = true }
      }, [tab, activeId, stepsKey, obsColumnsKey, analysesNonce, analysesPayload])

      // 嵌入图：tab 可见 + embed 步骤完成时才请求
      React.useEffect(function () {
        if (tab !== 'embedding' || !activeId) return undefined
        if (doneSteps.indexOf('embed') < 0) {
          setEmbedding(null)
          setEmbeddingError('')
          setEmbeddingLoading(false)
          return undefined
        }
        if (colorBy.kind === 'gene' && !String(colorBy.key || '').trim()) {
          setEmbedding(null)
          setEmbeddingError('')
          setEmbeddingLoading(false)
          return undefined
        }
        var cancelled = false
        var controller = typeof AbortController === 'function' ? new AbortController() : null
        var requestColor = colorBy
        setEmbeddingLoading(true)
        setEmbeddingError('')
        clientApi.getEmbedding(activeId, requestColor, controller ? { signal: controller.signal } : undefined)
          .then(function (payload) {
            if (cancelled) return
            var data = normalizeEmbedding(payload)
            if (!data.nPoints) {
              setEmbedding(null)
              setEmbeddingLoading(false)
              return
            }
            setEmbedding({ key: stepsKey + '|' + colorValueOf(requestColor) + '|' + embedNonce, data: data })
            setHighlight(null)
            setEmbeddingLoading(false)
          })
          .catch(function (error) {
            if (cancelled || isAbortError(error)) return
            setEmbedding(null)
            setEmbeddingError(formatError(error) || T('embedFailed'))
            setEmbeddingLoading(false)
          })
        return function () {
          cancelled = true
          if (controller) {
            try {
              controller.abort()
            } catch (error) {
              /* ignore */
            }
          }
        }
      }, [tab, activeId, stepsKey, colorBy.kind, colorBy.key, embedNonce])

      // 标记基因：标记/聚类 tab 打开且 cluster 步骤完成时按需拉（服务端有缓存）
      React.useEffect(function () {
        if ((tab !== 'markers' && tab !== 'clusters') || !activeId) return undefined
        if (doneSteps.indexOf('cluster') < 0) {
          setMarkersLoading(false)
          setMarkersError('')
          return undefined
        }
        var key = stepsKey + '|' + markersTop + '|' + markersMethod
        if (markers && markers.key === key) return undefined
        var cancelled = false
        setMarkersLoading(true)
        setMarkersError('')
        clientApi.markers(activeId, { top_n: markersTop, method: markersMethod }).then(function (payload) {
          if (cancelled) return
          setMarkers({ key: key, groups: normalizeMarkersGroups(payload) })
          setMarkersLoading(false)
        }).catch(function (error) {
          if (cancelled || isAbortError(error)) return
          setMarkersError(formatError(error) || T('requestFailed'))
          setMarkersLoading(false)
        })
        return function () { cancelled = true }
      }, [tab, activeId, stepsKey, markersTop, markersMethod, markersNonce, markers])

      // 聚类柱状图：聚类 tab 打开时自动生成一次（同 steps 复用）
      React.useEffect(function () {
        if (tab !== 'clusters' || !activeId || doneSteps.indexOf('cluster') < 0) return undefined
        if (clusterBar && clusterBar.key === stepsKey) return undefined
        var cancelled = false
        setClusterBarLoading(true)
        setClusterBarError('')
        clientApi.plot(activeId, 'cluster_bar', {}).then(function (payload) {
          if (cancelled) return
          var figure = normalizeFigure(payload)
          if (!figure) throw new Error(T('clusterBarFailed'))
          setClusterBar({ key: stepsKey, figure: figure })
          setClusterBarLoading(false)
        }).catch(function (error) {
          if (cancelled || isAbortError(error)) return
          setClusterBarError(formatError(error) || T('clusterBarFailed'))
          setClusterBarLoading(false)
        })
        return function () { cancelled = true }
      }, [tab, activeId, stepsKey, clusterBar, clusterBarNonce])

      // 日志 tab：打开时拉一次 /activity
      React.useEffect(function () {
        if (tab !== 'log' || !activeId) return undefined
        var cancelled = false
        setLogLoading(true)
        setLogError('')
        clientApi.getActivity().then(function (payload) {
          if (cancelled) return
          var entries = payload && Array.isArray(payload.entries)
            ? payload.entries
            : (payload && Array.isArray(payload.activity) ? payload.activity : [])
          setLog(normalizeActivity(entries))
          setLogLoading(false)
        }).catch(function (error) {
          if (cancelled || isAbortError(error)) return
          setLogError(formatError(error) || T('requestFailed'))
          setLogLoading(false)
        })
        return function () { cancelled = true }
      }, [tab, activeId, logNonce])

      // ---- 通用动作 ----

      function loadHealth() {
        if (healthRequestedRef.current) return Promise.resolve(null)
        healthRequestedRef.current = true
        return clientApi.getHealth().then(function (payload) {
          var next = normalizeEngine(payload && payload.engine ? payload.engine : payload)
          setHealthEngine(next)
          return next
        }).catch(function () {
          healthRequestedRef.current = false
          return null
        })
      }

      function refreshState(silent) {
        return clientApi.getState().then(function (payload) {
          var next = normalizeState(payload)
          setState(next)
          setBootError('')
          setPollError('')
          return next
        }).catch(function (error) {
          if (silent) setPollError(formatError(error) || T('pollFailed'))
          else setBootError(formatError(error) || T('requestFailed'))
          return null
        })
      }

      async function runTask(key, task) {
        setBusy(function (previous) {
          var next = Object.assign({}, previous)
          next[key] = true
          return next
        })
        setErrors(function (previous) {
          var next = Object.assign({}, previous)
          next[key] = ''
          return next
        })
        try {
          return await task()
        } catch (error) {
          if (!isAbortError(error)) {
            setErrors(function (previous) {
              var next = Object.assign({}, previous)
              next[key] = formatError(error) || T('requestFailed')
              return next
            })
          }
          return null
        } finally {
          setBusy(function (previous) {
            var next = Object.assign({}, previous)
            delete next[key]
            return next
          })
        }
      }

      function recordStepSummary(stepId, payload) {
        var summary = payload && payload.summary && typeof payload.summary === 'object' ? payload.summary : null
        if (!summary) return
        if (!summary.title && !(Array.isArray(summary.bullets) && summary.bullets.length)) return
        setStepSummaries(function (previous) {
          var next = Object.assign({}, previous)
          next[stepId] = summary
          return next
        })
      }

      function copyText(key, text) {
        var clipboard = null
        try {
          clipboard = typeof navigator !== 'undefined' && navigator ? navigator.clipboard : null
        } catch (error) {
          clipboard = null
        }
        if (!clipboard || typeof clipboard.writeText !== 'function') {
          setCopyError(T('copyFailed'))
          return
        }
        try {
          var promise = clipboard.writeText(String(text === null || text === undefined ? '' : text))
          if (promise && typeof promise.then === 'function') {
            promise.then(function () {
              setCopyError('')
              setCopiedKey(key)
              if (copyTimerRef.current) clearTimeout(copyTimerRef.current)
              copyTimerRef.current = setTimeout(function () {
                copyTimerRef.current = null
                setCopiedKey('')
              }, 1600)
            }).catch(function () {
              setCopyError(T('copyFailed'))
            })
          } else {
            setCopyError('')
            setCopiedKey(key)
          }
        } catch (error) {
          setCopyError(T('copyFailed'))
        }
      }

      // ---- 顶栏动作 ----

      async function handleDemo() {
        var payload = await runTask('demo', function () { return clientApi.openDemo() })
        if (payload) await refreshState(true)
      }

      async function handleOpenPath(filePath) {
        var payload = await runTask('open', function () { return clientApi.openPath(filePath) })
        if (payload) await refreshState(true)
      }

      async function handleUpload(file) {
        if (!file) return
        var name = typeof file.name === 'string' ? file.name : ''
        if (!/\.(h5ad|h5)$/i.test(name)) {
          setErrors(function (previous) {
            var next = Object.assign({}, previous)
            next.top = T('h5adNeeded')
            return next
          })
          return
        }
        var payload = await runTask('upload', function () { return clientApi.upload(file) })
        if (payload) await refreshState(true)
      }

      async function handleSelect(datasetId) {
        if (!datasetId || datasetId === activeId) return
        var payload = await runTask('select', function () { return clientApi.selectDataset(datasetId) })
        if (payload) await refreshState(true)
      }

      async function handleRefreshAll() {
        healthRequestedRef.current = false
        loadHealth()
        setEmbedNonce(function (value) { return value + 1 })
        setMarkersNonce(function (value) { return value + 1 })
        setClusterBarNonce(function (value) { return value + 1 })
        setLogNonce(function (value) { return value + 1 })
        await refreshState(false)
      }

      // ---- 流程动作 ----

      async function handleRunStep(stepId, params) {
        if (!activeId) return
        var payload = await runTask('step:' + stepId, function () {
          return clientApi.runStep(activeId, stepId, params)
        })
        if (!payload) return
        recordStepSummary(stepId, payload)
        setEmbedNonce(function (value) { return value + 1 })
        setMarkersNonce(function (value) { return value + 1 })
        setClusterBarNonce(function (value) { return value + 1 })
        await refreshState(true)
      }

      async function handleSaveAnnotation(mapping) {
        if (!activeId) return null
        var payload = await runTask('annotate', function () {
          return clientApi.runStep(activeId, 'annotate', { mapping: mapping })
        })
        if (!payload) return null
        recordStepSummary('annotate', payload)
        setAnnotationSaved(T('annotationSaved'))
        setAdopted(function (previous) {
          var next = Object.assign({}, previous)
          Object.keys(mapping).forEach(function (clusterId) { next[clusterId] = true })
          return next
        })
        setEmbedNonce(function (value) { return value + 1 })
        setMarkersNonce(function (value) { return value + 1 })
        await refreshState(true)
        return payload
      }

      function handleChangeAnnotation(clusterId, value) {
        setAnnotationSaved('')
        setAnnotation(function (previous) {
          var next = Object.assign({}, previous)
          next[clusterId] = value
          return next
        })
      }

      async function handleAdoptSuggestions(suggestionList) {
        var mapping = buildAnnotationMapping(suggestionList)
        if (!Object.keys(mapping).length) return
        var merged = Object.assign({}, annotation)
        Object.keys(mapping).forEach(function (clusterId) { merged[clusterId] = mapping[clusterId] })
        setAnnotation(merged)
        await handleSaveAnnotation(merged)
      }

      async function handleGenerateQc() {
        if (!activeId) return
        if (!qcMetrics.length) return
        var payload = await runTask('qc', function () {
          return clientApi.plot(activeId, 'qc_violin', { metrics: qcMetrics.slice() })
        })
        if (!payload) return
        var figure = normalizeFigure(payload)
        if (!figure) return
        setQcFigure(figure)
      }

      async function handleDotplot() {
        if (!activeId) return
        var payload = await runTask('dotplot', function () {
          return clientApi.plot(activeId, 'dotplot', { top_n: 5 })
        })
        if (!payload) return
        var figure = normalizeFigure(payload)
        if (figure) setDotplot(figure)
      }

      async function handleGenerateClusterBar() {
        if (!activeId) return
        var payload = await runTask('clusterBar', function () {
          return clientApi.plot(activeId, 'cluster_bar', {})
        })
        if (!payload) return
        var figure = normalizeFigure(payload)
        if (figure) setClusterBar({ key: stepsKey, figure: figure })
      }

      function handleReloadMarkers() {
        setMarkers(null)
        setMarkersNonce(function (value) { return value + 1 })
      }

      // ---- 记录 tab 动作（v2）----

      async function handleAudit() {
        if (!activeId) return
        var payload = await runTask('audit', function () { return clientApi.getAudit(activeId) })
        if (payload) setAudit(normalizeAudit(payload))
      }

      async function handleAutoannotate() {
        if (!activeId) return
        var payload = await runTask('autoann', function () { return clientApi.getAutoannotate(activeId) })
        if (payload) setAutoann(normalizeAutoannotate(payload))
      }

      async function handleScan() {
        if (!activeId) return
        var parsed = parseResolutions(scanInput)
        if (parsed.error) {
          setErrors(function (previous) {
            var next = Object.assign({}, previous)
            next.scan = parsed.error
            return next
          })
          return
        }
        var payload = await runTask('scan', function () {
          return clientApi.resolutionScan(activeId, parsed.list)
        })
        if (payload) setScan(normalizeResolutionScan(payload))
      }

      async function handleExport(what) {
        if (!activeId) return
        setExportWhat(what)
        var payload = await runTask('export', function () { return clientApi.exportData(activeId, what) })
        setExportWhat('')
        if (!payload) return
        setExportResult(normalizeExportResult(payload))
      }

      // ---- 预设流程 ----

      async function handleRunRecipe() {
        if (!activeId || !recipeId) return
        setRecipeResult(null)
        var payload = await runTask('recipe', function () {
          return clientApi.runRecipe(activeId, recipeId)
        })
        if (!payload) return
        setRecipeResult({
          recipe: typeof payload.recipe === 'string' ? payload.recipe : recipeId,
          title: typeof payload.title === 'string' ? payload.title : recipeId,
          completed: Number(payload.completed) || 0,
          total: Number(payload.total) || 0,
        })
        if (Array.isArray(payload.steps)) {
          setStepSummaries(function (previous) {
            var next = Object.assign({}, previous)
            payload.steps.forEach(function (item) {
              if (!item || typeof item !== 'object') return
              var stepId = String(item.step || '')
              if (!STEP_BY_ID[stepId]) return
              next[stepId] = {
                title: typeof item.title === 'string' ? item.title : '',
                bullets: Array.isArray(item.bullets) ? item.bullets : [],
              }
            })
            return next
          })
        }
        setEmbedNonce(function (value) { return value + 1 })
        setMarkersNonce(function (value) { return value + 1 })
        setClusterBarNonce(function (value) { return value + 1 })
        await refreshState(true)
      }

      // ---- 嵌入图动作 ----

      function handleColorChange(value) {
        if (value === 'gene') {
          setColorBy({ kind: 'gene', key: String(geneQuery || '').trim() })
          return
        }
        var parts = String(value || '').split(':')
        var kind = parts[0] === 'obs' ? 'obs' : 'cluster'
        var key = parts.length > 1 ? parts.slice(1).join(':') : (kind === 'cluster' ? 'clusters' : 'n_genes')
        setColorBy({ kind: kind, key: key })
      }

      function handleGeneSubmit() {
        var gene = String(geneQuery || '').trim()
        if (!gene) return
        setColorBy({ kind: 'gene', key: gene })
      }

      // ---- v3：导入动作 ----

      function openImportDialog() {
        setImportError('')
        setImportCandidates(EMPTY_LIST)
        setImportCandidatesDisabled(false)
        setInspectResult(null)
        setInspectError('')
        setDialogNonce(function (value) { return value + 1 })
        setDialogOpen(true)
      }

      function closeImportDialog() {
        setDialogOpen(false)
        setImportError('')
        setImportCandidates(EMPTY_LIST)
        setImportCandidatesDisabled(false)
        setInspectResult(null)
        setInspectError('')
      }

      function reloadFormats() {
        setFormats(null)
        setFormatsError('')
      }

      async function handleImport(options) {
        if (importBusy) return
        setImportBusy(true)
        setImportError('')
        setInspectError('')
        setImportCandidates(EMPTY_LIST)
        setImportCandidatesDisabled(false)
        lastImportRef.current = options
        try {
          var opened = null
          if (options.sourceType === 'file') {
            if (!options.file) throw new Error(T('importNeedSource'))
            opened = await clientApi.upload(options.file)
          } else {
            var importBody = buildImportBody(options)
            if (!importBody.path && !importBody.url && !(importBody.paths && importBody.paths.length)) {
              throw new Error(T('importNeedSource'))
            }
            opened = await clientApi.importDataset(options)
          }
          await finishOpen(opened)
        } catch (error) {
          if (error && error.code === 'ambiguous') {
            var formatted = formatCandidates(error)
            setImportCandidates(formatted.items)
            setImportCandidatesDisabled(options.sourceType === 'file')
            setImportError(formatError(error) || T('importCandidates'))
          } else {
            setImportError(formatError(error) || T('requestFailed'))
          }
        } finally {
          setImportBusy(false)
        }
      }

      function handlePickCandidate(candidate) {
        if (!candidate) return
        var last = lastImportRef.current
        if (!last || last.sourceType === 'file') return
        var patch = candidate.patch || {}
        var next = Object.assign({}, last)
        if (patch.member !== undefined) next.member = patch.member
        if (patch.object !== undefined) next.object = patch.object
        handleImport(next)
      }

      async function handleInspect(options) {
        if (inspectBusy) return
        setInspectBusy(true)
        setInspectError('')
        setInspectResult(null)
        try {
          var payload = await clientApi.inspectSource(options)
          setInspectResult(normalizeImportInspect(payload))
        } catch (error) {
          setInspectError(formatError(error) || T('requestFailed'))
        } finally {
          setInspectBusy(false)
        }
      }

      // ---- v3：高级分析动作 ----

      function handleSelectAnalysis(item) {
        if (!item || !item.id) return
        setSelectedAnalysis(item)
        setAnalysisValues(paramDefaults(item.params, obsColumns, groupColumn))
        setAnalysisResult(null)
        setAnalysisError('')
      }

      function handleAnalysisValueChange(key, value) {
        setAnalysisValues(function (previous) {
          var next = Object.assign({}, previous)
          next[key] = value
          return next
        })
      }

      async function handleRunAnalysis() {
        if (!activeId || !selectedAnalysis || selectedAnalysis.external) return
        var built = collectAnalysisParams(selectedAnalysis.params, analysisValues)
        if (built.error) {
          setAnalysisError(built.error)
          return
        }
        setAnalysisError('')
        setAnalysisResult(null)
        var payload = await runTask('analysis', function () {
          return clientApi.runAnalysis(activeId, selectedAnalysis.id, built.params)
        })
        if (!payload) return
        setAnalysisResult(normalizeAnalysisResult(payload))
        setPreviewNonce(function (value) { return value + 1 })
        await refreshState(true)
      }

      function handleReloadAnalyses() {
        analysesKeyRef.current = ''
        setAnalysesPayload(null)
        setAnalysesError('')
        setAnalysesNonce(function (value) { return value + 1 })
      }

      // ---- v4：样本 & 分组动作 ----

      function handleSampleName(sampleId, value) {
        setSampleDirty(true)
        setSampleRows(function (previous) {
          return previous.map(function (row) {
            return row.id === sampleId ? Object.assign({}, row, { name: value }) : row
          })
        })
      }

      function handleSampleGroup(sampleId, value) {
        setSampleDirty(true)
        setSampleRows(function (previous) {
          return previous.map(function (row) {
            return row.id === sampleId ? Object.assign({}, row, { group: value }) : row
          })
        })
      }

      function handleAddGroup() {
        var value = String(groupDraft || '').trim()
        if (!value) return
        setSampleGroups(function (previous) {
          return previous.indexOf(value) >= 0 ? previous : previous.concat([value]).sort()
        })
        setGroupDraft('')
      }

      function collectSampleGroups(rows, groups) {
        var list = normalizeStringList(groups)
        for (var index = 0; index < rows.length; index++) {
          var group = String(rows[index].group || '').trim()
          if (group && list.indexOf(group) < 0) list.push(group)
        }
        list.sort()
        return list
      }

      async function saveSamples(rows, groups) {
        if (!activeId || !Array.isArray(rows) || !rows.length) return null
        var body = buildSamplesUpdate(rows, {
          groups: groups,
          sampleColumn: obsPayload ? obsPayload.sampleColumn : 'sample',
          groupColumn: groupColumn,
        })
        var payload = await runTask('samplesSave', function () {
          return clientApi.samplesUpdate(activeId, body)
        })
        if (!payload) return null
        var nextSamples = normalizeSamples(payload.samples)
        if (nextSamples.length) {
          setSampleRows(parseSampleRows({ samples: nextSamples, groups: payload.groups }))
        }
        var nextGroups = normalizeStringList(payload.groups)
        if (nextGroups.length) setSampleGroups(nextGroups)
        setSampleDirty(false)
        var updatedCells = Number(payload.updated_cells)
        setSamplesMessage(
          T('samplesSaved') +
          (Number.isFinite(updatedCells) ? ' · ' + tf(T, 'samplesUpdatedCells', formatNumber(updatedCells)) : ''),
        )
        setPreviewNonce(function (value) { return value + 1 })
        setSamplesNonce(function (value) { return value + 1 })
        await refreshState(true)
        return payload
      }

      function handleSaveSamples() {
        return saveSamples(sampleRowsState, sampleGroups)
      }

      async function handleAutoGroupRun() {
        if (!activeId) return
        var pattern = String(autoPattern || '').trim()
        if (!pattern) return
        var payload = await runTask('samplesAuto', function () {
          return clientApi.samplesAutoGroup(activeId, pattern)
        })
        if (!payload) return
        setSampleMapPreview(null)
        setAutoPreview(parseAutoGroupMapping(payload))
      }

      async function handleAutoGroupApply() {
        if (!autoPreview || !autoPreview.rows.length) return
        var mapping = {}
        for (var index = 0; index < autoPreview.rows.length; index++) {
          mapping[autoPreview.rows[index].id] = autoPreview.rows[index].group
        }
        var nextRows = applySampleGroupMapping(sampleRowsState, mapping)
        var nextGroups = collectSampleGroups(nextRows, sampleGroups)
        setSampleRows(nextRows)
        setSampleGroups(nextGroups)
        setSampleDirty(true)
        setAutoPreview(null)
        await saveSamples(nextRows, nextGroups)
      }

      async function handleImportMapRun() {
        if (!activeId) return
        var filePath = String(sampleMapPath || '').trim()
        if (!filePath) return
        var payload = await runTask('samplesMap', function () {
          return clientApi.samplesImportMap(activeId, filePath)
        })
        if (!payload) return
        setAutoPreview(null)
        setSampleMapPreview(normalizeImportMapResult(payload))
      }

      async function handleImportMapApply() {
        if (!sampleMapPreview) return
        var nextRows = applySampleGroupMapping(sampleRowsState, sampleMapPreview.updates)
        var nextGroups = collectSampleGroups(nextRows, sampleGroups)
        setSampleRows(nextRows)
        setSampleGroups(nextGroups)
        setSampleDirty(true)
        setSampleMapPreview(null)
        await saveSamples(nextRows, nextGroups)
      }

      async function handleJoinObsRun() {
        if (!activeId) return
        var filePath = String(joinPath || '').trim()
        if (!filePath) return
        var payload = await runTask('obsJoin', function () {
          return clientApi.obsJoin(activeId, filePath, joinKey)
        })
        if (!payload) return
        setJoinResult(normalizeJoinObsResult(payload))
        setSamplesNonce(function (value) { return value + 1 })
        setPreviewNonce(function (value) { return value + 1 })
        await refreshState(true)
      }

      var samplesPanelNode = activeId
        ? h(SamplesPanel, {
            rows: sampleRows,
            groups: sampleGroups,
            groupOptions: sampleGroupOptions,
            groupColors: scanGroupColors,
            statusText: sampleStatusText,
            dirty: sampleDirty,
            loading: obsLoading,
            error: obsError,
            saving: Boolean(busy.samplesSave),
            notice: samplesMessage,
            groupDraft: groupDraft,
            onGroupDraft: setGroupDraft,
            onAddGroup: handleAddGroup,
            onNameChange: handleSampleName,
            onGroupChange: handleSampleGroup,
            onSave: handleSaveSamples,
            onDeleteGroup: handleDeleteSampleGroup,
            onRemoveSample: handleRemoveSample,
            onReload: function () { setSamplesNonce(function (value) { return value + 1 }) },
            autoPattern: autoPattern,
            autoPreview: autoPreview,
            autoBusy: Boolean(busy.samplesAuto),
            onAutoPattern: setAutoPattern,
            onAutoRun: handleAutoGroupRun,
            onAutoApply: handleAutoGroupApply,
            onAutoCancel: function () { setAutoPreview(null) },
            importPath: sampleMapPath,
            importPreview: sampleMapPreview,
            mapBusy: Boolean(busy.samplesMap),
            onImportPath: setSampleMapPath,
            onImportRun: handleImportMapRun,
            onImportApply: handleImportMapApply,
            onImportCancel: function () { setSampleMapPreview(null) },
            joinPath: joinPath,
            joinKey: joinKey,
            joinResult: joinResult,
            joinBusy: Boolean(busy.obsJoin),
            onJoinPath: setJoinPath,
            onJoinKey: setJoinKey,
            onJoinRun: handleJoinObsRun,
          })
        : null

      // ---- v5：扫描 / GEO / 确认导入 动作 ----

      function finishOpen(opened) {
        setDialogOpen(false)
        setPreviewNonce(function (value) { return value + 1 })
        setSamplesNonce(function (value) { return value + 1 })
        var openedDataset = opened && opened.dataset && typeof opened.dataset === 'object' ? opened.dataset : {}
        if (openedDataset.import_state) {
          setOpenedImportState({
            datasetId: openedDataset.id || '',
            state: openedDataset.import_state,
            notes: openedDataset.notes || EMPTY_LIST,
          })
        }
        var openedSamples = Array.isArray(openedDataset.samples) ? openedDataset.samples : []
        if (openedSamples.length > 1) {
          setMergeNotice({
            count: openedSamples.length,
            name: openedDataset.name || openedDataset.id || '',
          })
          setTab('record')
        } else {
          setTab('embedding')
        }
        return refreshState(true)
      }

      async function handlePickPath() {
        if (scanLoading) return
        setScanError('')
        try {
          var payload = await clientApi.pickPath('dir', T('scanBrowse'))
          if (!payload || payload.ok !== true) {
            setScanError((payload && payload.error) || T('scanBrowseUnavailable'))
            return
          }
          if (payload.cancelled) return
          var first = Array.isArray(payload.paths) ? String(payload.paths[0] || '') : ''
          if (first) setScanRoot(first)
          else setScanError(T('scanBrowseUnavailable'))
        } catch (error) {
          setScanError(formatError(error) || T('scanBrowseUnavailable'))
        }
      }

      async function handleScanRun() {
        if (scanLoading) return
        var root = String(scanRoot || '').trim()
        if (!root) {
          setScanError(T('scanRootRequired'))
          return
        }
        setScanLoading(true)
        setScanError('')
        setScanConfirmError('')
        try {
          var payload = await clientApi.scanDataRoot(root, scanRecursive !== false)
          var normalizedScan = normalizeScanResult(payload)
          setScanResult(normalizedScan)
          var nextRows = createScanRows(normalizedScan)
          setScanRows(function (previous) {
            var groupByKey = {}
            for (var index = 0; index < previous.length; index++) {
              var row = previous[index]
              if (row.path && row.group) groupByKey[row.path + '|' + row.trioPrefix] = row.group
            }
            return nextRows.map(function (row) {
              var key = row.path + '|' + row.trioPrefix
              return Object.assign({}, row, { group: groupByKey[key] || '' })
            })
          })
          if (normalizedScan.errors.length) setScanError(normalizedScan.errors.join('；'))
        } catch (error) {
          setScanError(formatError(error) || T('requestFailed'))
        } finally {
          setScanLoading(false)
        }
      }

      function updateScanRow(key, patch) {
        setScanRows(function (previous) {
          return previous.map(function (row) {
            return row.key === key ? Object.assign({}, row, patch) : row
          })
        })
      }

      function handleScanRowSelect(key, selected) {
        updateScanRow(key, { selected: selected === true })
      }

      function handleScanRowName(key, name) {
        updateScanRow(key, { name: value })
      }

      function handleScanRowGroup(key, group) {
        updateScanRow(key, { group: group })
      }

      function handleScanSelectAll() {
        setScanRows(function (previous) { return setScanRowsSelected(previous, true) })
      }

      function handleScanClearAll() {
        setScanRows(function (previous) { return setScanRowsSelected(previous, false) })
      }

      function handleAddScanGroup() {
        var name = String(scanNewGroup || '').trim()
        if (!name) return
        setScanGroups(function (previous) {
          return previous.indexOf(name) >= 0 ? previous : previous.concat([name])
        })
        setScanGroupColors(function (previous) {
          return Object.assign({}, previous, { [name]: scanNewGroupColor || '#4e79a7' })
        })
        setScanNewGroup('')
      }

      function handleApplyBatchGroup() {
        var group = String(scanBatchGroup || '').trim()
        if (!group) return
        setScanRows(function (previous) { return applyGroupToRows(previous, group, true) })
      }

      async function handleGeoQuery() {
        var gse = String(geoGse || '').trim()
        if (!gse || geoLoading) return
        setGeoLoading(true)
        setGeoError('')
        try {
          var payload = await clientApi.queryGeo(gse)
          var geo = normalizeGeoResult(payload)
          var rows = Array.isArray(scanRows) ? scanRows : EMPTY_LIST
          var matched = matchGeoToSamples(rows, payload)
          if (rows.length) setScanRows(matched.rows)
          setGeoResult({
            gse: geo.gse || gse.toUpperCase(),
            title: geo.title,
            samples: geo.samples,
            unmatched: matched.unmatched,
            matchedCount: matched.matched,
          })
        } catch (error) {
          setGeoError(formatError(error) || T('requestFailed'))
        } finally {
          setGeoLoading(false)
        }
      }

      async function handleConfirmSamples() {
        var validation = validateSampleSelection(scanRows)
        if (!validation.ok) {
          setScanConfirmError(validation.errors.join('；'))
          return
        }
        setScanConfirmError('')
        setImportBusy(true)
        try {
          var payload = buildSamplesOpenPayload(scanRows)
          var opened = await clientApi.importDataset(payload)
          await finishOpen(opened)
        } catch (error) {
          setScanConfirmError(formatError(error) || T('requestFailed'))
        } finally {
          setImportBusy(false)
        }
      }

      // ---- v6：降维聚类动作 ----

      function currentDimredSettings() {
        return normalizeDimredSettings({
          engine: dimEngine,
          harmony: dimHarmony,
          batchKey: dimBatchKey,
          resolution: dimResolution,
          nPcs: dimNpcs,
          neighbors: dimNeighbors,
          nTopGenes: dimTopGenes,
        })
      }

      function autoBatchKey() {
        var explicit = String(dimBatchKey || '').trim()
        if (explicit) return explicit
        for (var index = 0; index < obsColumns.length; index++) {
          if (obsColumns[index].key === 'sample') return 'sample'
        }
        return ''
      }

      function applyDimredRunResult(payloads) {
        var bullets = []
        for (var index = 0; index < payloads.length; index++) {
          var payload = payloads[index]
          var summary = payload && payload.summary && typeof payload.summary === 'object' ? payload.summary : null
          if (summary && Array.isArray(summary.bullets)) bullets = bullets.concat(summary.bullets)
        }
        setDimWarnings(pickDimredWarnings(activeDataset ? activeDataset.warnings : EMPTY_LIST, [{ summary: { bullets: bullets } }]))
        setEmbedNonce(function (value) { return value + 1 })
        setMarkersNonce(function (value) { return value + 1 })
        setPreviewNonce(function (value) { return value + 1 })
      }

      async function handleDimredRecipe() {
        if (!activeId) return
        var payload = await runTask('dimredRecipe', function () {
          return clientApi.runRecipe(activeId, DIMRED_RECIPE_ID)
        })
        if (!payload) return
        var payloads = []
        var bullets = []
        var steps = Array.isArray(payload.steps) ? payload.steps : []
        for (var index = 0; index < steps.length; index++) {
          var item = steps[index]
          if (!item || typeof item !== 'object') continue
          if (Array.isArray(item.bullets)) bullets = bullets.concat(item.bullets)
          payloads.push({ summary: { bullets: Array.isArray(item.bullets) ? item.bullets : [] } })
        }
        setStepSummaries(function (previous) {
          var next = Object.assign({}, previous)
          for (var stepIndex = 0; stepIndex < steps.length; stepIndex++) {
            var stepItem = steps[stepIndex]
            if (!stepItem || !STEP_BY_ID[String(stepItem.step || '')]) continue
            next[stepItem.step] = {
              title: stepItem.title || stepTitleOf(stepItem.step),
              bullets: Array.isArray(stepItem.bullets) ? stepItem.bullets : [],
            }
          }
          return next
        })
        setDimWarnings(pickDimredWarnings(activeDataset ? activeDataset.warnings : EMPTY_LIST, payloads))
        setEmbedNonce(function (value) { return value + 1 })
        setMarkersNonce(function (value) { return value + 1 })
        setPreviewNonce(function (value) { return value + 1 })
        await refreshState(true)
      }

      async function handleDimredRecluster() {
        if (!activeId) return
        var settings = currentDimredSettings()
        var batchKey = autoBatchKey()
        var result = await runTask('dimredRecluster', async function () {
          var first = await clientApi.runStep(activeId, 'cluster', buildReclusterParams(settings, batchKey))
          var second = await clientApi.runStep(activeId, 'embed', buildReclusterEmbedParams(settings, batchKey))
          return { first: first, second: second }
        })
        if (!result) return
        recordStepSummary('cluster', result.first)
        recordStepSummary('embed', result.second)
        applyDimredRunResult([result.first, result.second])
        await refreshState(true)
      }

      // ---- v6：手动注释动作 ----

      function handleSelectCluster(clusterId) {
        setSelectedCluster(String(clusterId))
        setSelectedGene('')
        setAnnotPlotFigure(null)
      }

      function handleSelectGene(gene) {
        setSelectedGene(String(gene || ''))
      }

      async function handlePlotMarker(kind) {
        if (!activeId || !selectedGene) return
        var payload = await runTask('annotPlot', function () {
          return clientApi.plot(activeId, kind === 'heatmap' ? 'heatmap' : 'dotplot', {
            groupby: 'clusters',
            genes: [selectedGene],
            top_n: 5,
          })
        })
        if (!payload) return
        var figure = normalizeFigure(payload)
        if (figure) setAnnotPlotFigure(figure)
      }

      function handleFeaturePlot() {
        if (!selectedGene) return
        setColorBy({ kind: 'gene', key: selectedGene })
        setEmbedNonce(function (value) { return value + 1 })
        setTab('embedding')
      }

      function handleApplyAnnotationCsv() {
        var parsed = parseAnnotationCsv(annotCsvText)
        if (!parsed.rows.length) {
          setErrors(function (previous) {
            var next = Object.assign({}, previous)
            next.annotate = T('annotCsvEmpty')
            return next
          })
          return
        }
        setAnnotation(function (previous) { return Object.assign({}, previous, parsed.mapping) })
        setAnnotationSaved(tf(T, 'annotCsvApplied', parsed.rows.length))
        setAnnotCsvOpen(false)
        setAnnotCsvText('')
      }

      // ---- v6：分组维护动作 ----

      async function handleDeleteSampleGroup(group) {
        if (isGroupUsedBySamples(sampleRowsState, group)) {
          setSamplesMessage(T('samplesGroupUsed'))
          return
        }
        var nextGroups = removeGroupFromList(sampleGroups, group)
        setSampleGroups(nextGroups)
        var payload = await saveSamples(sampleRowsState, nextGroups)
        if (payload) setSamplesMessage(T('samplesGroupDeleted'))
      }

      async function handleRemoveSample(sampleId) {
        if (!activeId) return
        var payload = await runTask('samplesRemove', function () {
          return clientApi.samplesRemove(activeId, [sampleId])
        })
        if (!payload) return
        setSamplesMessage(tf(T, 'samplesRemoved', payload.removed_cells || 0))
        await refreshState(true)
      }

      // ---- v6：Agent 意图动作 ----

      function findAnalysisById(id) {
        var list = analysesPayload ? analysesPayload.analyses : EMPTY_LIST
        for (var index = 0; index < list.length; index++) {
          if (list[index].id === id) return list[index]
        }
        return null
      }

      function firstMarkerGene() {
        if (!markerGroups.length) return ''
        var genes = markerGroups[0].genes || []
        return genes.length ? String(genes[0].gene || '') : ''
      }

      function handleAgentSubmit() {
        var text = String(agentInput || '').trim()
        if (!text) return
        var parsed = parseAgentIntent(text, {
          steps: doneSteps,
          hasMarkers: markerGroups.length > 0,
          hasAnnotation: Object.keys(annotationMapping).length > 0,
          hasGroups: sampleGroups.length > 0,
        })
        setAgentSuggestions(parsed.suggestions)
        setAgentFallback(parsed.fallback)
        setAgentResults({})
      }

      function handleAgentCopy(suggestion) {
        if (!suggestion) return
        copyText('agent:' + suggestion.id, agentActionInstruction(suggestion.action))
      }

      async function runAgentAction(action) {
        if (!activeId || !action) throw new Error('没有可执行动作')
        var type = String(action.type || '')
        if (type === 'recipe') {
          var recipePayload = await clientApi.runRecipe(activeId, action.id || DIMRED_RECIPE_ID)
          var recipeSteps = recipePayload && Array.isArray(recipePayload.steps) ? recipePayload.steps : []
          return '预设流程完成 ' + recipeSteps.length + ' 步'
        }
        if (type === 'recluster') {
          var settings = normalizeDimredSettings(Object.assign({}, currentDimredSettings(), action.params || {}))
          var batchKey = autoBatchKey()
          var clusterPayload = await clientApi.runStep(activeId, 'cluster', buildReclusterParams(settings, batchKey))
          var embedPayload = await clientApi.runStep(activeId, 'embed', buildReclusterEmbedParams(settings, batchKey))
          recordStepSummary('cluster', clusterPayload)
          recordStepSummary('embed', embedPayload)
          var clusterSummary = clusterPayload && clusterPayload.summary && Array.isArray(clusterPayload.summary.bullets)
            ? clusterPayload.summary.bullets
            : []
          return clusterSummary.length ? clusterSummary.join('；') : ('已重聚类（resolution=' + settings.resolution + '）')
        }
        if (type === 'markers') {
          var markersPayload = await clientApi.markers(activeId, { top_n: markersTop, method: markersMethod })
          setMarkers({ key: stepsKey + '|' + markersTop + '|' + markersMethod, groups: normalizeMarkersGroups(markersPayload) })
          setTab('markers')
          return '标记基因已更新'
        }
        if (type === 'plot') {
          var plotPayload = await clientApi.plot(activeId, action.kind === 'heatmap' ? 'heatmap' : 'dotplot', {
            groupby: 'clusters',
            genes: action.gene ? [action.gene] : [],
            top_n: 5,
          })
          var plotFigure = normalizeFigure(plotPayload)
          if (plotFigure) setDotplot(plotFigure)
          setTab('markers')
          return '已生成 ' + String(action.kind || 'dotplot')
        }
        if (type === 'featureplot') {
          var gene = String(action.gene || selectedGene || firstMarkerGene() || '')
          if (!gene) throw new Error('没有可用基因：先在标记基因 tab 选一个基因')
          setColorBy({ kind: 'gene', key: gene })
          setEmbedNonce(function (value) { return value + 1 })
          setTab('embedding')
          return '已按 ' + gene + ' 给嵌入图着色'
        }
        if (type === 'annotate') {
          var mapping = Object.keys(annotationMapping).length ? annotationMapping : {}
          if (!Object.keys(mapping).length) throw new Error('还没有填写细胞类型注释')
          await handleSaveAnnotation(mapping)
          return '注释已保存'
        }
        if (type === 'analysis') {
          var item = findAnalysisById(action.id)
          var params = item ? paramDefaults(item.params, obsColumns, groupColumn) : {}
          var analysisPayload = await clientApi.runAnalysis(activeId, action.id, params)
          var analysisNormalized = normalizeAnalysisResult(analysisPayload)
          setAnalysisResult(analysisNormalized)
          if (item) setSelectedAnalysis(item)
          setTab('analysis')
          return analysisNormalized.summary.bullets[0] || '分析完成'
        }
        if (type === 'audit') {
          var auditPayload = await clientApi.getAudit(activeId)
          var auditData = normalizeAudit(auditPayload)
          setAudit(auditData)
          setTab('record')
          return '自检得分 ' + auditData.score
        }
        if (type === 'export') {
          var exportPayload = await clientApi.exportData(activeId, action.what || 'report')
          var exported = normalizeExportResult(exportPayload)
          setExportResult(exported)
          setTab('record')
          return '已导出：' + exported.file
        }
        throw new Error('不支持的动作：' + type)
      }

      async function handleAgentRun(suggestion) {
        if (!suggestion || agentBusyId) return
        var id = suggestion.id
        setAgentBusyId(id)
        setAgentResults(function (previous) {
          var next = Object.assign({}, previous)
          next[id] = ''
          return next
        })
        try {
          var summary = await runAgentAction(suggestion.action)
          setAgentResults(function (previous) {
            var next = Object.assign({}, previous)
            next[id] = summary || T('agentDone')
            return next
          })
          setEmbedNonce(function (value) { return value + 1 })
          setPreviewNonce(function (value) { return value + 1 })
          await refreshState(true)
        } catch (error) {
          setAgentResults(function (previous) {
            var next = Object.assign({}, previous)
            next[id] = T('agentFailed') + '：' + (formatError(error) || T('requestFailed'))
            return next
          })
        } finally {
          setAgentBusyId('')
        }
      }

      // ---- 渲染 ----

      var tabNode = null
      if (tab === 'embedding') {
        tabNode = h(EmbeddingTab, {
          colorBy: colorBy,
          geneQuery: geneQuery,
          genes: preview ? preview.genes : EMPTY_LIST,
          layout: embedding && embedding.data ? embedding.data.layout : '',
          hasEmbed: doneSteps.indexOf('embed') >= 0,
          pcaDone: doneSteps.indexOf('pca') >= 0,
          running: runningStepIds.indexOf('embed') >= 0,
          loading: embeddingLoading,
          error: embeddingError,
          data: embedding ? embedding.data : null,
          highlight: highlight,
          onHighlight: setHighlight,
          onColorChange: handleColorChange,
          onGeneQuery: setGeneQuery,
          onGeneSubmit: handleGeneSubmit,
          colorOptions: colorOptions,
          onRefresh: function () { setEmbedNonce(function (value) { return value + 1 }) },
          onRunEmbed: function () { handleRunStep('embed', defaultParams('embed')) },
        })
      } else if (tab === 'qc') {
        tabNode = h(QcTab, {
          activeId: activeId,
          metrics: qcMetrics,
          loading: Boolean(busy.qc),
          error: errors.qc || '',
          figure: qcFigure,
          qcDone: doneSteps.indexOf('qc') >= 0,
          onToggleMetric: function (metricId) {
            setQcMetrics(function (previous) {
              var index = previous.indexOf(metricId)
              if (index >= 0) return previous.filter(function (item) { return item !== metricId })
              return previous.concat([metricId])
            })
          },
          onGenerate: handleGenerateQc,
        })
      } else if (tab === 'markers') {
        tabNode = h(MarkersTab, {
          activeId: activeId,
          clusterDone: doneSteps.indexOf('cluster') >= 0,
          groups: markerGroups,
          loading: markersLoading,
          error: markersError,
          topN: markersTop,
          method: markersMethod,
          onTopNChange: setMarkersTop,
          onMethodChange: setMarkersMethod,
          onReload: handleReloadMarkers,
          dotplot: dotplot,
          dotplotLoading: Boolean(busy.dotplot),
          dotplotError: errors.dotplot || '',
          onDotplot: handleDotplot,
          annotation: annotation,
          selectedCluster: selectedCluster,
          selectedGene: selectedGene,
          onSelectCluster: handleSelectCluster,
          onSelectGene: handleSelectGene,
          onChangeAnnotation: handleChangeAnnotation,
          onSaveAnnotation: function () { handleSaveAnnotation(annotationMapping) },
          annotationSaved: annotationSaved,
          annotationError: errors.annotate || '',
          savingAnnotation: Boolean(busy.annotate),
          plotBusy: Boolean(busy.annotPlot),
          plotError: errors.annotPlot || '',
          plotFigure: annotPlotFigure,
          onPlotMarker: handlePlotMarker,
          onFeaturePlot: handleFeaturePlot,
          csvOpen: annotCsvOpen,
          csvText: annotCsvText,
          onToggleCsv: function () { setAnnotCsvOpen(function (value) { return !value }) },
          onCsvText: setAnnotCsvText,
          onApplyCsv: handleApplyAnnotationCsv,
        })
      } else if (tab === 'clusters') {
        tabNode = h(ClustersTab, {
          activeId: activeId,
          clusterDone: doneSteps.indexOf('cluster') >= 0,
          groups: markerGroups,
          markersLoading: markersLoading,
          markersError: markersError,
          annotation: annotation,
          onChangeAnnotation: handleChangeAnnotation,
          hasMapping: annotationCount > 0,
          saving: Boolean(busy.annotate),
          saved: annotationSaved,
          saveError: errors.annotate || '',
          onSave: function () { handleSaveAnnotation(annotationMapping) },
          onReloadMarkers: handleReloadMarkers,
          barFigure: clusterBar ? clusterBar.figure : null,
          barLoading: clusterBarLoading || Boolean(busy.clusterBar),
          barError: clusterBarError || errors.clusterBar || '',
          onGenerateBar: handleGenerateClusterBar,
          dimredCard: h(DimredCard, {
            activeId: activeId,
            dim: dimredSettings,
            batchOptions: batchColumnOptions,
            onEngine: setDimEngine,
            onHarmony: setDimHarmony,
            onBatchKey: setDimBatchKey,
            onResolution: setDimResolution,
            onNPcs: setDimNpcs,
            onNeighbors: setDimNeighbors,
            onTopGenes: setDimTopGenes,
            runningRecipe: Boolean(busy.dimredRecipe),
            runningRecluster: Boolean(busy.dimredRecluster),
            warnings: dimWarnings,
            error: errors.dimredRecipe || errors.dimredRecluster || '',
            onRunRecipe: handleDimredRecipe,
            onRecluster: handleDimredRecluster,
          }),
        })
      } else if (tab === 'record') {
        tabNode = h(RecordTab, {
          activeId: activeId,
          samplesPanel: samplesPanelNode,
          history: history,
          clusterDone: doneSteps.indexOf('cluster') >= 0,
          audit: audit,
          auditLoading: Boolean(busy.audit),
          auditError: errors.audit || auditError,
          onAudit: handleAudit,
          autoann: autoann,
          autoannLoading: Boolean(busy.autoann),
          autoannError: errors.autoann || autoannError,
          onAutoannotate: handleAutoannotate,
          adopted: adopted,
          adopting: Boolean(busy.annotate),
          adoptError: errors.annotate || '',
          onAdoptCluster: function (suggestion) { handleAdoptSuggestions([suggestion]) },
          onAdoptAll: function () { handleAdoptSuggestions(autoann ? autoann.suggestions : EMPTY_LIST) },
          scan: scan,
          scanLoading: Boolean(busy.scan),
          scanError: errors.scan || scanError,
          scanInput: scanInput,
          onScanInputChange: setScanInput,
          onScan: handleScan,
          exportResult: exportResult,
          exportLoadingKind: busy.export ? exportWhat : '',
          exportError: errors.export || exportError,
          onExport: handleExport,
          copiedKey: copiedKey,
          onCopyPath: function () { copyText('export-path', exportResult ? exportResult.file : '') },
        })
      } else if (tab === 'log') {
        tabNode = h(LogTab, {
          entries: log || activity,
          loading: logLoading,
          error: logError,
          onRefresh: function () { setLogNonce(function (value) { return value + 1 }) },
        })
      } else {
        tabNode = h(AnalysisTab, {
          activeId: activeId,
          groups: analysesGroups,
          loading: analysesLoading,
          listError: analysesError,
          selected: selectedAnalysis,
          doneSteps: doneSteps,
          obsColumns: obsColumns,
          values: analysisValues,
          result: analysisResult,
          error: analysisError || errors.analysis || '',
          running: Boolean(busy.analysis),
          onSelect: handleSelectAnalysis,
          onReload: handleReloadAnalyses,
          onValueChange: handleAnalysisValueChange,
          onRun: handleRunAnalysis,
        })
      }

      var recipeStrip = null
      var recipeTitle = ''
      if (recipeProgress) {
        recipeTitle = recipeProgress.title || recipeProgress.recipe || ''
      } else if (recipeId) {
        var matchedRecipe = recipes.filter(function (recipe) { return recipe.id === recipeId })
        recipeTitle = matchedRecipe.length ? matchedRecipe[0].title : recipeId
      }
      var recipeText = ''
      if (busy.recipe) {
        recipeText = (recipeTitle ? recipeTitle + ' · ' : '') +
          (recipeProgress ? recipeProgressText(recipeProgress) : T('recipeRunning'))
      } else if (recipeResult) {
        recipeText = tf(T, 'recipeDone', recipeResult.completed) + (recipeResult.title ? ' · ' + recipeResult.title : '')
      }
      if (recipeText) {
        recipeStrip = h('div', { className: 'sc-recipe-strip' },
          busy.recipe ? h(Spinner) : null,
          h('span', null, recipeText))
      }

      var topActionError = errors.top || errors.open || errors.demo || errors.upload || errors.select || errors.recipe || ''
      var centerNode = null
      if (booting) {
        centerNode = LoadingLine({ text: T('loading') })
      } else if (!activeId) {
        centerNode = h(EmptyState, { busy: busy, onDemo: handleDemo, onOpenImport: openImportDialog })
      } else {
        centerNode = tabNode
      }

      return h('div', { className: 'sc-root' },
        h('style', { 'data-plugin': 'dsh-singlecell' }, CSS),
        h(PluginErrorBoundary, null,
          h(TopBar, {
            datasets: normalized.datasets,
            active: activeId,
            engine: engine,
            recipes: recipes,
            recipeId: recipeId,
            recipeRunning: Boolean(busy.recipe),
            onRecipeChange: setRecipeId,
            onRunRecipe: handleRunRecipe,
            busy: busy,
            onSelect: handleSelect,
            onOpenImport: openImportDialog,
            onDemo: handleDemo,
            onRefresh: handleRefreshAll,
          })),
        bootError
          ? h('div', { className: 'sc-top-err' },
              h(ErrorBanner, {
                message: bootError,
                onRetry: function () { refreshState(false) },
                retryLabel: T('retry'),
              }))
          : null,
        topActionError
          ? h('div', { className: 'sc-top-err' }, h('div', { className: 'sc-inline-err' }, topActionError))
          : null,
        pollError
          ? h('div', { className: 'sc-top-err' }, h('div', { className: 'sc-warn-text' }, T('pollFailed') + '：' + pollError))
          : null,
        recipeStrip,
        mergeNotice
          ? h('div', { className: 'sc-merge-notice' },
              h('span', { className: 'sc-merge-notice-text' },
                tf(T, 'importMergedPrefix', formatNumber(mergeNotice.count)) +
                (mergeNotice.name ? '：' + mergeNotice.name : '') +
                ' · ' + T('importMergedHint')),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm sc-btn-primary',
                onClick: function () { setMergeNotice(null); setTab('record') },
              }, T('importGoSamples')),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm',
                onClick: function () { setMergeNotice(null) },
              }, '×'))
          : null,
        importStateVisible
          ? h('div', { className: 'sc-merge-notice' },
              h('span', { className: 'sc-merge-notice-text' }, importStateHintText),
              h('button', {
                type: 'button',
                className: 'sc-btn sc-btn-sm',
                onClick: function () {
                  setImportStateDismissed(function (previous) {
                    return Object.assign({}, previous, { [activeId || '']: true })
                  })
                },
              }, '×'))
          : null,
        h('div', { className: 'sc-body' },
          h(PluginErrorBoundary, null,
            h(PipelinePanel, {
              activeId: activeId,
              doneSteps: doneSteps,
              runningStepIds: runningStepIds,
              summaries: mergedSummaries,
              errors: errors,
              annotationMapping: annotationMapping,
              annotationCount: annotationCount,
              onRunStep: handleRunStep,
            })),
          h(PluginErrorBoundary, null,
            h('div', { className: 'sc-center' },
              h(TabsBar, { tab: tab, onChange: setTab }),
              h('div', { className: 'sc-pane' }, centerNode))),
          h(PluginErrorBoundary, null,
            h(AgentPanel, {
              engine: engine,
              activeDataset: activeDataset,
              suggestions: suggestions,
              activity: activity,
              doneSteps: doneSteps,
              copiedKey: copiedKey,
              copyError: copyError,
              onCopy: copyText,
              agentInput: agentInput,
              onAgentInput: setAgentInput,
              onAgentSubmit: handleAgentSubmit,
              agentSuggestions: agentSuggestions,
              agentResults: agentResults,
              agentBusyId: agentBusyId,
              agentFallback: agentFallback,
              onAgentRun: handleAgentRun,
              onAgentCopy: handleAgentCopy,
            }))),
        dialogOpen
          ? h(ImportDialog, {
              key: 'import-' + dialogNonce,
              scanRoot: scanRoot,
              scanRecursive: scanRecursive,
              scanLoading: scanLoading,
              scanError: scanError,
              scanResult: scanResult,
              scanRows: scanRows,
              scanGroups: scanGroups,
              scanGroupColors: scanGroupColors,
              scanNewGroup: scanNewGroup,
              scanNewGroupColor: scanNewGroupColor,
              scanBatchGroup: scanBatchGroup,
              scanConfirmError: scanConfirmError,
              geoGse: geoGse,
              geoLoading: geoLoading,
              geoResult: geoResult,
              geoError: geoError,
              onScanRootChange: setScanRoot,
              onPickPath: handlePickPath,
              onScanRecursive: setScanRecursive,
              onScan: handleScanRun,
              onScanRowSelect: handleScanRowSelect,
              onScanRowName: handleScanRowName,
              onScanRowGroup: handleScanRowGroup,
              onSelectAll: handleScanSelectAll,
              onClearAll: handleScanClearAll,
              onScanNewGroupChange: setScanNewGroup,
              onScanNewGroupColorChange: setScanNewGroupColor,
              onAddScanGroup: handleAddScanGroup,
              onBatchGroupChange: setScanBatchGroup,
              onApplyBatchGroup: handleApplyBatchGroup,
              onGeoGseChange: setGeoGse,
              onGeoQuery: handleGeoQuery,
              onConfirmSamples: handleConfirmSamples,
              onScanError: setScanError,
              formats: formats,
              formatsError: formatsError,
              busy: importBusy,
              inspectBusy: inspectBusy,
              error: importError,
              inspectError: inspectError,
              candidates: importCandidates,
              candidatesDisabled: importCandidatesDisabled,
              inspectResult: inspectResult,
              onImport: handleImport,
              onInspect: handleInspect,
              onPickCandidate: handlePickCandidate,
              onLocalError: setImportError,
              onClose: closeImportDialog,
              onReloadFormats: reloadFormats,
            })
          : null)
    }

    function WorkbenchPage() {
      return h(PluginErrorBoundary, null, h(WorkbenchBody))
    }

    // ==========================================================================
    // 插件入口：locale 字典 + 两个槽
    // ==========================================================================

    function apply(ctx) {
      var locale = ctx.locale
      if (locale && typeof locale.register === 'function') {
        var registerDictionaries = function () {
          try {
            return locale.register(NS, { zh: DICT_ZH, en: DICT_EN })
          } catch (error) {
            warnPlugin('failed to register locale dictionaries', error)
            return undefined
          }
        }
        try {
          if (typeof ctx.effect === 'function') ctx.effect(registerDictionaries, 'dsh-singlecell: locale')
          else registerDictionaries()
        } catch (error) {
          warnPlugin('failed to install locale effect', error)
        }
        T = createTranslator(locale)
      }

      var slots = ctx.slots
      if (!slots && typeof ctx.get === 'function') {
        try {
          slots = ctx.get('slots')
        } catch (error) {
          slots = null
        }
      }
      if (!slots || typeof slots.register !== 'function' || typeof slots.inject !== 'function') {
        warnPlugin('slots service unavailable; UI not registered')
        return
      }

      // inject/register 的返回值可能是函数，也可能是 undefined：统一收拢，卸载时逐个清掉。
      var disposers = []
      function keep(handle) {
        if (typeof handle === 'function') disposers.push(handle)
        return handle
      }

      try {
        if (typeof ctx.effect === 'function') {
          ctx.effect(function () {
            return function () {
              while (disposers.length) {
                var dispose = disposers.pop()
                try {
                  dispose()
                } catch (error) {
                  /* 单个 disposer 失败不影响其它 */
                }
              }
            }
          }, 'dsh-singlecell: slot disposers')
        }
      } catch (error) {
        warnPlugin('failed to install slot cleanup', error)
      }

      try {
        keep(slots.inject('main', function () {
          return keep(slots.register({ name: 'main', key: 'singlecell', locale: NS }, WorkbenchPage))
        }))
      } catch (error) {
        warnPlugin('failed to register main slot', error)
      }

      try {
        keep(slots.inject('sidebar.panellist', function () {
          return keep(slots.register({
            name: 'sidebar.panellist',
            id: 'singlecell',
            order: 3,
            label: function () { return T('nav') },
            locale: NS,
          }, SidebarIcon))
        }))
      } catch (error) {
        warnPlugin('failed to register sidebar.panellist slot', error)
      }
    }

    // __test__ 只是自测/联调口（宿主不消费）：纯函数 + 组件引用。
    return {
      name: NS,
      inject: ['slots', 'locale'],
      apply: apply,
      __test__: Object.assign({}, SC_TEST, {
        components: {
          WorkbenchPage: WorkbenchPage,
          WorkbenchBody: WorkbenchBody,
          TopBar: TopBar,
          EmptyState: EmptyState,
          TabsBar: TabsBar,
          PipelinePanel: PipelinePanel,
          StepCard: StepCard,
          AgentPanel: AgentPanel,
          EmbeddingTab: EmbeddingTab,
          EmbeddingCanvas: EmbeddingCanvas,
          QcTab: QcTab,
          MarkersTab: MarkersTab,
          ClustersTab: ClustersTab,
          RecordTab: RecordTab,
          LogTab: LogTab,
          ImportDialog: ImportDialog,
          AnalysisTab: AnalysisTab,
          SamplesPanel: SamplesPanel,
          DimredCard: DimredCard,
          AnnotationWorkbench: AnnotationWorkbench,
          EngineBadge: EngineBadge,
          SidebarIcon: SidebarIcon,
          PluginErrorBoundary: PluginErrorBoundary,
        },
      }),
    }
  }

  // ==========================================================================
  // 顶层守卫：node 里只暴露纯函数给单测；浏览器里注册 Client 模块。
  // ==========================================================================

  if (typeof window === 'undefined') {
    globalThis.__SC_TEST__ = SC_TEST
  } else if (window.__ModuleLoader__ && typeof window.__ModuleLoader__.load === 'function') {
    window.__ModuleLoader__.load({ id: NS, factory: factory })
  }
})()
