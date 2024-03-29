// bdfファイルの読み込み
Merge "fem.bdf";

// メッシュの生成アルゴリズムの指定
Mesh.Algorithm = 6; // Frontalアルゴリズムを使用

// メッシュオプションの設定
//Mesh.ElementType = 2; // 三角形要素を指定


// メッシュの生成
Mesh 2;

// あるプロパティIDに対応する要素のグループを作成
//PhysicalGroupID = AddPhysicalGroup(2, {要素IDのリスト});
//SetPhysicalName(2, PhysicalGroupID, "プロパティ名");

//PhysicalGroupID = AddPhysicalGroup(2, {3});
//SetPhysicalName(2, PhysicalGroupID, "aluminumplate");

//PhysicalGroupID2 = AddPhysicalGroup(2, {4, 5, 6, 7, 8, 9, 10, 11, 12, 13,15,16,17,18,19});
//SetPhysicalName(2, PhysicalGroupID2, "rib");

// メッシュ形式の指定（バージョン2のMSHフォーマット）
Mesh.MshFileVersion = 2.2;
// メッシュのエクスポート
Save "converted_mesh2.msh";