# UNLSI
## イントロダクション
全てMATLABによって記述された非構造メッシュ低次Source-Doubletパネル法のソルバーです。 
構造用のソルバーとしてshell-FEM(準3次元薄板要素有限要素法）を実装しているため、空力解析の結果を用いて静的な構造解析や空力弾性解析も行うことができます。
その他航空機の機体設計に有用なプログラムを目指して様々な機能が実装されています。
また、このUNLSIを用いて機体形状の最適化を行えるUNGRADEも付属しています。

このソフトウェアの使用方法については、2024年のUNLSIワークショップで用いた資料を含めていますのでそちらを参考してください。

V7.01よりMATLABのApp Designerを用いたGUIアプリケーションを追加しました。こちらはMATLABがなくても（誰でも）Runtimeをインストールすることで使用することができます。
GUI版はReleaseページで添付している.exeを用いてインストールしてください。（現在はWindowsのみ）

UNLSIGUIフォルダには、アプリケーションをコンパイルする前の.mlapp（とGUI用のUNLSI.m、他バージョンとの競合を防ぐため）が入っています。

## 特徴(UNLSI)
### パネル法による空力解析
- MATLABのクラスによって記述されているため使いやすく、また継承によって他のプログラムへの組み込みが簡単。
- openVSP (https://github.com/OpenVSP/OpenVSP) の機体記述形式であるvspgeomに直接対応しているため、機体形状作成が簡単。
- openVSPの機能の一つであるCFDmesh機能によって作成した非構造メッシュにて解析が行えるため、メッシュの状態が良く、解析精度が良い。
- 差分法による動安定微係数の推算が可能。舵効きも簡易的ではあるが推算できる。
- Actuator Diskによってプロペラ後流を含めた解析ができる。
### パネル法と数値微分法による航空機運動解析
- パネル法の特徴を利用した高速な空力安定微係数推算が実行できる。
- 空力係数の補間曲面を作成することでシミュレータへのデータの受け渡しが簡単にできる。
- 各運動モードの固有振動数や減衰係数も計算できる。
### Shell-FEMによる構造解析
- openVSPのFEAmeshやGMSH等で作成できる構造用のメッシュに対し、薄肉要素の解析ができる。
- 空力構造メッシュ間の変換にはメッシュ変形手法を用いて、各メッシュにおける物理量を写像している。
- 質量行列と減衰行列を用いた空力弾性解析ができる。またモード解析にも対応しており、モード空間での空力弾性解析もできる。
  
## 特徴(UNGRADE)
- openVSPのdesファイルに最適化したい変数を登録するだけで機体形状が最適化される。
- massPropatiesなどの各種解析にも対応。
- openVSPではなく自作の機体形状作成プログラム等でも動作する。
- 機体形状の更新に様々なアルゴリズムを実装済みのため、問題に合わせて調整できる。
- UNLSIの動安定微係数解析や舵効きの解析にも対応。固有振動数やつり合い迎え角等を評価関数・制約条件に含めて最適化できる。
- 構造解析の連成も可能。かけられる計算コストに応じて弱連成（基準形状の空力荷重で構造解析を行う）、強連成（構造解析後の構造変形を考慮して空力解析を行い、釣り合い点を求める）の選択ができる。

## インストール
1. 適当なところにcloneもしくはダウンロードする。
2. MATLABにて当該フォルダ(UNLSI)にパスを通す。（UNLSIのみでよい）

The MIT License
Copyright (c) 2023 Naoto Morita

以下に定める条件に従い、本ソフトウェアおよび関連文書のファイル（以下「ソフトウェア」）の複製を取得するすべての人に対し、ソフトウェアを無制限に扱うことを無償で許可します。これには、ソフトウェアの複製を使用、複写、変更、結合、掲載、頒布、サブライセンス、および/または販売する権利、およびソフトウェアを提供する相手に同じことを許可する権利も無制限に含まれます。

上記の著作権表示および本許諾表示を、ソフトウェアのすべての複製または重要な部分に記載するものとします。

ソフトウェアは「現状のまま」で、明示であるか暗黙であるかを問わず、何らの保証もなく提供されます。ここでいう保証とは、商品性、特定の目的への適合性、および権利非侵害についての保証も含みますが、それに限定されるものではありません。 作者または著作権者は、契約行為、不法行為、またはそれ以外であろうと、ソフトウェアに起因または関連し、あるいはソフトウェアの使用またはその他の扱いによって生じる一切の請求、損害、その他の義務について何らの責任も負わないものとします。


Copyright 2023 Naoto Morita

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

また、本ソフトウェアに付帯するソフトウェアとしてOpenVSPおよびdistanceVertex2Mesh.mを用いています。
上記ソフトウェアのライセンス表記を以下に示します。

openVSP

Copyright (c) 2012 United States Government as represented by the Administrator for The National Aeronautics and Space Administration. All Rights Reserved.


DISTANCEVERTEX2MESH - calculate the distance between vertices and a mesh

Author: Christopher Haccius

Telecommunications Lab, Saarland University, Germany

email: haccius@nt.uni-saarland.de

March 2015; Last revision: 26-March-2015

