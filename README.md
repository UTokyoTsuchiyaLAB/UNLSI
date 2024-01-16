# UNLSI
## イントロダクション
全てMATLABによって記述された非構造メッシュ低次Source-Doubletパネル法のソルバーです。 
また、このUNLSIを用いて機体形状の最適化を行えるUNGRADEも付属しています。
## 特徴(UNLSI)
- MATLABのクラスによって記述されているため使いやすく、また継承によって他のプログラムへの組み込みが簡単。
- openVSP (https://github.com/OpenVSP/OpenVSP) の機体記述形式であるvspgeomに直接対応しているため、機体形状作成が簡単。
- openVSPの機能の一つであるCFDmesh機能によって作成した非構造メッシュにて解析が行えるため、メッシュの状態が良く、解析精度が良い。
- 差分法による動安定微係数の推算が可能。舵効きも簡易的ではあるが推算可能。
- Actuator Diskによってプロペラ後流を含めた解析が可能。
- v5.02よりshellFEM機能を実装した。openVSPのFEAmeshで作成したコンポーネントのメッシュに対し、薄肉要素の解析ができる。
- 各メッシュ間の変換にはメッシュ変形手法を用いて、各メッシュにおける物理量を写像している。
## 特徴(UNGRADE)
- ※under development
- openVSPのdesファイルに最適化したい変数を登録するだけで機体形状が最適化される。
- massPropatiesなどの各種解析にも対応。
- openVSPではなく自作の機体形状作成プログラム等でも動作する。
- 機体形状の更新に様々なアルゴリズムを実装済みのため、問題に合わせて調整できる。
- UNLSIの動安定微係数解析や舵効きの解析にも対応。固有振動数やつり合い迎え角等を評価関数・制約条件に含めて最適化できる。

## インストール
1. 適当なところにcloneもしくはダウンロードする。
2. MATLABにて当該フォルダ(UNLSI)にパスを通す。（UNLSIのみでよい）

## パブリックメソッド
### UNLSI
#### UNLSI(verts,connectivity,surfID,wakelineID,halfmesh)
            %%%%%%%%%%%Constructor%%%%%%%%%%%%%
            %verts:頂点座標
            %conectivity:各パネルの頂点ID
            %surfID:パネルのID
            %wakelineID:wakeをつけるエッジ上の頂点ID配列を要素にもつセル配列
            %halfmesh:半裁メッシュの場合に1を指定　省略した場合はy座標の値をみて自動判定するが、halfmeshでも-0.000があったりするので注意。
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### setREFS(obj,SREF,BREF,CREF)
            %%%%%%%%%%%Referentials Setting%%%%%%%%%%%%%
            %SREF:基準面積
            %BREF:横方向基準長(eg.翼幅)
            %CREF:縦方向基準長(eg.MAC,機体全長
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### setRotationCenter(obj,XYZREF)
            %%%%%%%%%%%Rotation Center Setting%%%%%%%%%%%%%
            %回転中心を設定する。
            %回転中心：モーメント計算の基準位置,主流角速度の回転中心
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### setCf(obj,flowNo,Re)
            %%%%%%%%%%%%%%%%%%%%%Cf計算における設定%%%%%%%%%%%%%%%%%%%
            %Cf計算の詳細はhttps://openvsp.org/wiki/doku.php?id=parasitedragが詳しい
            %flowNo:流れのID
            %Re:レイノルズ数
            %Lch:代表長さ
            %k:skin roughness value
                %Camouflage paint on aluminium
                %k = 1.015*(10^-5); 
                %Smooth paint
                %k = 0.634*(10^-5);
                %Produciton sheet metal
                %k = 0.405*(10^-5);
                %Polished sheet metal
                %k = 0.152*(10^-5);
                %Smooth molded composite
                %k = 0.052*(10^-5);
            %LTratio:層流の割合
            %coefficient:調整用の係数(デフォルト:1)

            %TODO:パネルIDによってCfの値を切り替えられるようにする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### solveFlow(obj,alpha,beta,Mach,Re)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %結果はobj.AERODATAに格納される。
            % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt 7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14:E(翼効率) 15:CFx 16:CFy 17:CFz 18:CMx 19:CMy 20:CMz 21:CMl 22:CMm 23:CMn 24:FOpt 
            %上記で求めていないものは0が代入される
            %flowNo:解きたい流れのID
            %alpha:迎角[deg]
            %beta:横滑り角[deg]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### solvePertPotential(obj,flowNo,alpha,beta,omega,propCalcFlag)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %Adjoint法実装のためポテンシャルの変動値のみ求める。

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### solveFlowForAdjoint(obj,u,flowNo,alpha,beta,omega,propCalcFlag)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %ポテンシャルから力を求める
            %結果は配列に出力される
            % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CDo 7:CDi 8:CDtot 9:CDt 10:CDtot_t 11:CS 12:L/D E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
            %上記で求めていないものは0が代入される
            %u: doubletの強さ
            %flowNo:解きたい流れのID
            %alpha:迎角[deg]
            %beta:横滑り角[deg]
            %omega:主流の回転角速度(deg/s)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### setDeflAngle(obj,ID,rotAxis,angle)
            %%%%%%%%%%%Rotation Center Setting%%%%%%%%%%%%%
            %疑似的な舵角を設定する。
            %そのID上のパネルの法線ベクトルをロドリゲスの回転ベクトルによって回転する
            % 誤差については要検討。
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### plotGeometry(obj,figureNo,triColor,colorlim,method,extrapmethod)
            %%%%%%%%%%%plotting%%%%%%%%%%%%%
            %機体メッシュもしくは状態量をプロットする。
            %カラー値は各パネルごとの値（例えばCp）を入れると、sccatteredInterpを用いて接点の値に補間しプロットする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            


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

