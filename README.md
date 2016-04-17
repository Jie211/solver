Iterative Krylov-Subspace Solvers
====
*-last update 16/4/17-*

**行列変換ツールのPythonバージョンを追加**

**日本語での説明を追加しました**
## News
* 2016/4/17 Add MatrixTransformTools, rewrite by Python2.
* 2016/4/12 Now Working on GMRES Method.
* 2016/3/20 Now Working on CUDA version Kskip-CG.

## Branches
* **cuda branche** is for parallel with CUDA.(now work at here)
* **makefile branche** is for make this solvers as a linux package by using GNUAutomake.(**last update 2015/12**)
* **cmdline branche** is a dev branche.(now not use)
* **master branche** is **not using** now.

## Description
Iterative Krylov-Subspace Solvers for large systems of linear equations. 
All the code is write be pure C.
This repository is part of my project.

## About Kskip method
NOLTA2015(English) [keynote](https://www.dropbox.com/s/ni0gt1m93izdhem/NOLTA2015_12_3.key?dl=0)

## About Variable Preconditioned method
Abe, Kuniyoshi, and Shao-Liang Zhang. "A variable preconditioning using the SOR method for GCR-like methods." Int. J. Numer. Anal. Model 2.2 (2005): 147-161.

Ikuno, S., Kawaguchi, Y., Fujita, N., Itoh, T., Nakata, S., & Watanabe, K. (2012). Iterative solver for linear system obtained by edge element: variable preconditioned method with mixed precision on GPU. Magnetics, IEEE Transactions on, 48(2), 467-470.

## Now Support method
* Conjugate Gradient method(CG).
* Conjugate Residual method(CR).
* Generalized Conjugate Residual method(GCR).
* Generalized Minimal REsidual method(GMRES).(Anyway, implementation done, meybe rewrite itlater.)
* Kskip method
	- Kskip-CG.(implementation done, but algorithm is not stability)
	- Kskip-CR.(implementation done, stil working...)
* Variable Preconditioned method(a.k.a VP)
	- VP Conjugate Gradient method(VPCG).
	- VP Conjugate Residual method(VPCR).
	- VP Generalized Conjugate Residual method(VPGCR).
  - VP Generalized Minimal REsidual method(VPGMRES).(Anyway, implementation done, meybe rewrite itlater.)

## Requirement
Hard

* CPU(of course...)
* NVIDIA GPU(capability version >= 2.x, architecture >= Kepler like NVIDIA GTX Titan)

Soft

* GCC compiler in CentOS need >= 4.4.7.(see CUDA Toolkit Documentation 1.1 System Requirements)
* CUDA Ver>6.0.
* OpenMP.

## Licence
Apache License

see ./LICENSE

Copyright (c) 2015-2016 GongChen <g211501428@edu.teu.ac.jp>

## ブランチ
* **cuda ブランチ** は今メインとして進んでいるブランチ、主に並列計算の部分をOpenMPからCUDAへ変換して作成しています。
* **makefile ブランチ** はGNUのAutomakeを利用して、configureとmake && make installでこの線形ソルバーをlinuxのプログラムとしてインストールきるように作成した。今は更新していない、新しい解法を全部完成する前にはまだ更新しないつもりです。(**last update 2015/12**)
* **cmdline ブランチ** は開発途中でのブランチ、今は使っていない.(now not use)
* **master ブランチ** は開発途中でのブランチ、今は使っていない.


## 説明
このリポジトリは連立一次方程式の**反復解法**である**クリロフ部分空間解法**のソルバー群である。

いわゆる**線形ソルバー**でです。

ソースコードは開発環境のために全部はCで書いています。

このリポッジトリは私の「並列化を前提とした通信回避Krylov部分空間解法」研究の一環となります。

## 通信回避解法について
こちら「非線形理論とその応用シンポジウム2015」のキーノートに参考してください　[keynote](https://www.dropbox.com/s/ni0gt1m93izdhem/NOLTA2015_12_3.key?dl=0)

## 可変的前処理付き解法について
下記の論文に参考してください

Abe, Kuniyoshi, and Shao-Liang Zhang. "A variable preconditioning using the SOR method for GCR-like methods." Int. J. Numer. Anal. Model 2.2 (2005): 147-161.

Ikuno, S., Kawaguchi, Y., Fujita, N., Itoh, T., Nakata, S., & Watanabe, K. (2012). Iterative solver for linear system obtained by edge element: variable preconditioned method with mixed precision on GPU. Magnetics, IEEE Transactions on, 48(2), 467-470.


## 実装した解法
* 共役勾配法　Conjugate Gradient method(CG).
* 共役残差法　Conjugate Residual method(CR).
* 一般共役残差法　Generalized Conjugate Residual method(GCR).
* 一般最小残差法　Generalized Minimal REsidual method(GMRES).(とりあえず実装した、これから書き直す可能性が高い.)
* 通信回避解法群　Kskip method
	- 通信回避共役勾配法　Kskip-CG.(時々収束が不安定)
	- 通信回避共役残差法　Kskip-CR.(アルゴリズムが再構成中)
* 可変的前処理付き解法群　Variable Preconditioned method(a.k.a VP)
	- 可変的前処理付き共役勾配法　VP Conjugate Gradient method(VPCG).
	- 可変的前処理付き共役残差法　VP Conjugate Residual method(VPCR).
	- 可変的前処理付き一般共役残差法　VP Generalized Conjugate Residual method(VPGCR).
  - 可変的前処理付き一般最小残差法　VP Generalized Minimal REsidual method(VPGMRES).(とりあえず実装した、これから書き直す可能性が高い.)

## 使い方
- ハードウェア編

	* マルチコアのCPU
	* NVIDIAのGPUが必要となる。
	* CUDAの開発環境がインストール済みのcapabilityが2.0以上、アーキテクチャがKepler世代或いはKepler世代以上が必要となる。

- ソフトウェア編

	* CUDAの開発環境としてGCCが必要となる、バージョンはOSによります、詳しい情報はCUDA Toolkit Documentation 1.1 System Requirementsに参照してください。

	* CUDAのバージョンが6.0より新しいものが必要となる.
	* OpenMPが必要となる.
- コンパイル
	* makefileを利用してコンパイルする(今後はAutomakeを作る予定)。

- サンプルの実行
	* 問題である行列を用意する。（下の二つの方法がある）
		- [ここ](https://www.dropbox.com/sh/cuvspwrca345www/AABVYFIuHmAAFQTH2vqkztmEa?dl=0)で最小テスト行列を使用する。
		- [MatrixMarket](http://math.nist.gov/MatrixMarket/)で行列をダウンロードする(解法により使う行列の種類が違う、CG CR kskipCG kskipCRに対しては正定値対称行列が必要、それ以外の解法では対称性行列が必要である。)
	
	* ダウンロードしたディレクトリ次のように配置する。
	
	~~~~
	./
	|-- Matrix
	|   `-- CSR
	|       `-- bcsstk14
	|-- solver
	|   |-- CRS
	|   |   |-- cg.cu
	|   |   |-- cg.h
	|   |   |-- cr.cu
	|   |   |-- cr.h
	|   |   |-- gcr.cu
	~~~~
	
	* オプション一覧は
	
	~~~~
	./a.out
	~~~~
	だけを実行する。

## ライセンス
Apache License

see ./LICENSE

Copyright (c) 2015-2016 GongChen <g211501428@edu.teu.ac.jp>

## Author
[Jie211@github](https://github.com/Jie211)

[blog](https://www.jie211.me)