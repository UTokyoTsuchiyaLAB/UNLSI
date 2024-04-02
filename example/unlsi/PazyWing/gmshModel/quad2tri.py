import meshio
import numpy as np

def quad_to_tri(quads):
    """四角形を二つの三角形に分割する"""
    tris = []
    for quad in quads:
        # 四角形[ABCD]を三角形[ABC]と[ACD]に分割
        tris.append([quad[0], quad[1], quad[2]])
        tris.append([quad[0], quad[2], quad[3]])
    return np.array(tris)

def process_mesh(input_filename, output_filename):
    # .mshファイルの読み込み
    mesh = meshio.read(input_filename)
    
    # 新しいセルのリストを作成
    new_cells = []
    for cell in mesh.cells:
        if cell.type == "quad":
            # 四角形セルを三角形に分割
            tris = quad_to_tri(cell.data)
            new_cells.append(meshio.CellBlock("triangle", tris))
        else:
            # 四角形以外のセルはそのまま追加
            new_cells.append(cell)
    
    # 新しいメッシュデータを作成
    new_mesh = meshio.Mesh(points=mesh.points, cells=new_cells)
    
    # 変換後のメッシュをファイルに保存
    meshio.write(output_filename, new_mesh, file_format="msh2")


# スクリプトの実行
input_filename = 'airfoiltest0110.msh'  # 元の.mshファイル名
output_filename = 'converted_mesh.msh'  # 出力される.mshファイル名
process_mesh(input_filename, output_filename)
