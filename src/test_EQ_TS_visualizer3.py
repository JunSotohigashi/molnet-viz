import EQ_TS_visualizer3
import rdkit.Chem

EQ_path = "manu_EQ_list.log"
TS_path = "manu_TS_list.log"

eq_block = """# Geometry of EQ 0, SYMMETRY = C2v 
C 	  -0.000112000000	   1.198714000000	   0.000000000000
H 	  -0.000007000000	   1.750081000000	   0.933062000000
H 	  -0.000007000000	   1.750081000000	  -0.933062000000
C 	   0.000000000000	  -0.182179000000	   0.000000000000
N 	   0.000098000000	  -1.371339000000	   0.000000000000
Energy    = -132.055704962587 (-132.055704962587 :    0.000000000000)
Spin(**2) =    0.770506670235
ZPVE      =    0.031274079270
Normal mode eigenvalues : nmode = 9
  0.005453116   0.006931664   0.019747598   0.042820693   0.043936614
  0.083153588   0.165421172   0.384847516   0.410750967"""


ts_block = """# Geometry of TS 0, SYMMETRY = Cs  
C 	   0.094288736605	   1.009283955327	   0.000012543986
H 	   0.176386821627	   1.541012952839	   0.937857062184
H 	   0.176296734553	   1.540896961790	  -0.937900534375
C 	  -0.842665550797	  -0.461183876690	  -0.000067286555
N 	   0.395665258015	  -0.484651993263	   0.000098214761
Energy    = -131.943386984554 (-131.943386984554 :    0.000000000000)
Spin(**2) =    0.758826152177
ZPVE      =    0.029651071075
Normal mode eigenvalues : nmode = 9
 -0.007783973   0.008953976   0.024170005   0.028146367   0.043815099
  0.084517047   0.119834834   0.384922643   0.419166250
CONNECTION : 0 - DC"""


def test_eq_load():
    eq1 = EQ_TS_visualizer3.EQ(eq_block)
    assert eq1.name == "EQ0"
    assert eq1.symmetry == "C2v"
    assert eq1.n_atoms == 5
    assert eq1.spin == 0.770506670235
    assert eq1.zpve == 0.031274079270
    assert eq1.nmode == 9
    assert eq1.eigenvalues == [
        0.005453116,
        0.006931664,
        0.019747598,
        0.042820693,
        0.043936614,
        0.083153588,
        0.165421172,
        0.384847516,
        0.410750967,
    ]
    

def test_bond():
    eq1 = EQ_TS_visualizer3.EQ(eq_block)
    eq1.determine_bonds()
    eq1_smiles = eq1.get_smiles()
    assert eq1_smiles == "[H]C([H])=C#N"
    



def test_ts_load():
    eq1 = EQ_TS_visualizer3.EQ(eq_block)
    eq_list = [eq1]
    ts1 = EQ_TS_visualizer3.TS(ts_block, eq_list)
    assert ts1.con_from == eq1
    assert ts1.con_to is None


def test_dc_judge():
    eq1 = EQ_TS_visualizer3.EQ(eq_block)
    eq_list = [eq1]
    ts1 = EQ_TS_visualizer3.TS(ts_block, eq_list)
    assert ts1.has_dc_connection()


test_bond()
