# 3DReconstruction
一些测试代码和日志之类的东西

### 2017/7/18
#### 思路
* 预处理（平滑Wiener filter）（需要实现）现在用的递归高斯平滑
* 图像分割（k-means clustering已尝试, iterative method, regional growth method, GASA见论文）（转换为二值图像）
* 形态学运算（用于去除二值化图像中的小空洞）（basic dilate）radius?
* 边缘检测 （cannyedge）
* 4-邻域标记（用matlab实现了）
* 肺轮廓补偿（跳过了）
* 提取肺实质（弄出来不对，只有一点点轮廓，占比8000/512/512）已解决：黑白反转（连通区域还是不对，左右以及胸腔全部都标记为同一个数了）已解决：改成4-邻域标记
* 提取ROI

#### Todo:
* 把matlab部分写成c++
* 把序列全部自动分割，然后三维重建
* ROI分割

### 2017/10/19
#### 已完成
* Input: a folder containing CT images (test folder contains 40 slices)
* Pipeline: input -> RecursiveGaussianImageFilter(~0.049s per slice) -> ScalarImageKmeansImageFilter(~4.0s per slice) -> CannyEdgeDetectionImageFilter(~0.45s per slice) -> MyConnectedRegionDetection(~0.09s per slice) -> itktovtk ->MarchingCubes
* Output: lung reconstruction result
#### Todo
* Segment skin and lung nodules
* improve efficiency
* test on a larger scale of data(~500 slices)

### 2017/10/22
#### 已完成
* Input: a folder containing CT images (test folder contains 40 slices)
* Pipeline: input -> BinaryThresholdImageFilter -> GrayscaleFillholeImageFilter -> BinaryMorphologicalOpeningImageFilter & BinaryErodeImageFilter -> SubtractImageFilter -> Segment from raw dicom file -> MarchingCubes
* Output: Skin reconstruction result
#### Todo
* improve efficiency of the algorithm for lung nodules segmentation
* test on a larger scale of data
* 最核心的挑战在哪里，创新点在哪里。别人的工作做到什么程度了。
* 在一张CT上点一个点或者画一个圈重建肿瘤区域。关键：生长规则的确定。
* 动态CT。每个窗拍摄1~2个呼吸周期，空间拼接。以后只需要采一次CT，根据结构光给出的表面数据就能对应肿瘤位置。
* 模糊图像情况下

2017/11/09
#### Todo
* 另一个病例跑一遍
* 表皮和肿瘤空洞问题优化
* 调研：找有没有相关工作 动态扫描
CT fluorescent(CT透视动态扫描)
IEEE transaction on medical image
Medical Image Analysis

elsevier science
ieee
springer

2017/12/16
#### 已完成
* Cine Scan病例yaoxin的头文件解析
* Skin segmentation
#### Todo
* 建模表面的运动模式 f(t;fixed_x)
* 运动最显著的区域
* 此病例没有肿瘤，另一个病例肿瘤分割，建模运动
#### 问题
* 数据有误，重新采集中

2017/12/30
#### 已完成
* 水平集方法分割结节
> GeodesicActiveContourLevelSetImageFilter接受两个输入，Input为水平集(LevelSetImage)，FeatureImage为根据图像梯度构造的图像。
* FeatureImage Pipeline: reader -> CurvatureAnisotropicDiffusionImageFilter -> GradientMagnitudeRecursiveGaussianImageFilter -> SigmoidImageFilter
* LevelSetImage Pipeline: reader -> fastMarchingImageFilter
* parameters: 
> seed(s) position, shortest distance from seed(s) to the nodule edge. 402 347 4.5 387 327 1.4
> sigmoid filter requires two parameters: $\alpha$ and $\beta$ -0.1 50.0
#### Todo
* 根据临近(z)分割出的结果计算后续分割的参数
* nodule segmentation in a series

2018/1/7
#### 已完成
* Level Set方法分割结节。输入initialSliceIndex, seedPosition, initialDistance可以从邻近Slice中分割结节
> 参数：
静态1 145 402 347 4.5
静态2 162 315 356 6.6
静态shenmingmou 176 374 303 8.4（GG）
* 问题：只依赖前一次结果更新seed风险太大了。非孤立型结节分不出来。