import meshio

# ファイルのパスを指定
file_path = r"C:\Users\aimit\git\AePW3-LDWG\00_Models\02_Pazy_TU_Delft\01_Models_GFEM\01_FlexibleTip\fem.bdf"

# 新しく生成されたファイルの保存先パス
output_path = r"C:\Users\aimit\git\AePW3-LDWG\00_Models\02_Pazy_TU_Delft\01_Models_GFEM\01_FlexibleTip\fem.msh"

# メッシュデータを読み込む
mesh = meshio.read(file_path)

# メッシュデータを書き出す
meshio.write(output_path, mesh)