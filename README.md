# README

本仓库`report.pdf`为课程设计报告，`code/`为所有文中运行的Matlab代码的全部。文件结构如下：

- fft: 所有设计的fft代码片段
  - fftNewPadded.m: Radix-2 CT算法的递归实现，含补零
  - fft_interative.m: Radix-2 CT算法的迭代实现
  - fft_radix4.m: Radix-4 CT算法的递归实现
  - fft_radix4_iterative.m: Radix-4 CT算法的迭代实现
  - stockham_fft.m: Stockham FFT算法的实现
  - fftdivpad.m: 分裂基FFT算法的实现
  - hybrid_fft.m: 混合串行：Hybrid_Serial FFT算法实现
  - fast_conv.m: 自行编写的基于FFT的快速卷积
- test: 本文中所有测试代码片段
  - test.m
  - hybridtest.m: 混合串行算法的代码测试
  - JudgeHybrid.m: 混合串行算法中，确认最佳N的代码
  - test_fast_conv.m: 快速卷积的测试


由于时间实在有限，以及分裂基算法编写的不够理想，所以最后我没得到预期的结果，也没得到理想的HS-FFT算法。因为从理论上分析应当是分裂基是最优的，但是实际实验效果和理论差别很大，我反复找代码哪里写错了，也没找到原因。真的特别遗憾，之前一直期望完成的并行优化也没能完成。希望等截止过后能继续做修改代码和实验，得到理想的实验结果。