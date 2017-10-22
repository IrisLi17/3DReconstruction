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


