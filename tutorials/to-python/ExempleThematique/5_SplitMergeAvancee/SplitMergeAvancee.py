# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Splitting and Merging a MED File Using the Advanced MEDLoader API
#
# ## Objective
#
# This exercise presents a complete and advanced use case of the advanced MEDLoader API. The goal is to create a multi-type mesh from scratch with 2 fields:
#
# - A cell field "CellField"
# - A node field "NodeField"
#
# We will then split these fields into two parts (for parallel processing by a code, for example) and also demonstrate how to merge two fields from disjoint pieces.
#
# ## Implementation Start
#
# Create an unstructured mesh `m0` from a structured mesh (meshDim=2, spaceDim=2) of 30*30. Each of the even cells of the mesh will be simplexized (i.e., split into triangles - using `MEDCouplingUMesh.simplexize(0)` method).

# +
import medcoupling as mc

m0 = mc.MEDCouplingCMesh()
arr = mc.DataArrayDouble(31, 1)
arr.iota(0.0)
m0.setCoords(arr, arr)
m0 = m0.buildUnstructured()
m00 = m0[::2]  # Extract even cells
m00.simplexize(0)
m01 = m0[1::2]
m0 = mc.MEDCouplingUMesh.MergeUMeshes([m00, m01])
m0.getCoords()[:] *= 1 / 15.0  # Illustrate how to quickly rescale a mesh
m0.setName("mesh")
# -

# <div class="alert alert-block alert-success">
# <b>Note:</b> setName() on "m0" is mandatory. Remember that in the MED file context, proper naming of meshes is fundamental.</div>
#
# Create the fields `cellField` and `nodeField` at timestep (5,6) and at time 5.6 s.

# Cell field
cellField = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
cellField.setTime(5.6, 5, 6)
cellField.setMesh(m0)
cellField.setName("CellField")
cellField.fillFromAnalytic(1, "exp(-((x-1)*(x-1)+(y-1)*(y-1)))")
cellField.getArray().setInfoOnComponent(0, "powercell [W]")
# Node field
nodeField = mc.MEDCouplingFieldDouble(mc.ON_NODES, mc.ONE_TIME)
nodeField.setTime(5.6, 5, 6)
nodeField.setMesh(m0)
nodeField.setName("NodeField")
nodeField.fillFromAnalytic(1, "exp(-((x-1)*(x-1)+(y-1)*(y-1)))")
nodeField.getArray().setInfoOnComponent(0, "powernode [W]")

# For example, this is obtained for "CellField":
#
# <img src="SplitAndMergeCell1.jpg" style="width:500px;">
#
# ## Mesh Partitioning
#
# Split `m0` into two distinct parts. The two parts will be named `proc0` and `proc1`. `proc0` will be the part inside the bounding box (`MEDCouplingUMesh.getCellsInBoundingBox()`) `[(0.,0.4),(0.,0.4)]` with a precision of 1e-10. `proc1` will be the complement (`DataArrayInt.buildComplement()`).

proc0 = m0.getCellsInBoundingBox([(0.0, 0.4), (0.0, 0.4)], 1e-10)
proc1 = proc0.buildComplement(m0.getNumberOfCells())

# <img src="SplitAndMergeCell2.jpg" style="width:500px;">
#
# ## Writing to 2 Separate MED Files
#
# Starting from the partitioning `proc0` and `proc1`, create 2 MED files named "proc0.med" and "proc1.med":

# +
nodeField0 = nodeField[proc0]
cellField0 = cellField[proc0]
cellField0.setMesh(nodeField0.getMesh())
nodeField1 = nodeField[proc1]
cellField1 = cellField[proc1]
cellField1.setMesh(nodeField1.getMesh())

proc0_fname = "proc0.med"
mc.WriteField(proc0_fname, nodeField0, True)
mc.WriteFieldUsingAlreadyWrittenMesh(proc0_fname, cellField0)

proc1_fname = "proc1.med"
mc.WriteField(proc1_fname, nodeField1, True)
mc.WriteFieldUsingAlreadyWrittenMesh(proc1_fname, cellField1)
# -

# Reading and merging the 2 separate MED files (less optimal)
#
# Starting from "proc0.med" and "proc1.med", read their respective "CellField" using the basic API, aggregate both and put the result in `cellField_read`:

cellField0_read = mc.ReadFieldCell("proc0.med", "mesh", 0, "CellField", 5, 6)
cellField1_read = mc.ReadFieldCell("proc1.med", "mesh", 0, "CellField", 5, 6)
cellField_read = mc.MEDCouplingFieldDouble.MergeFields(
    [cellField0_read, cellField1_read]
)

# <div class="alert alert-block alert-success">
# <b>Note:</b> It may seem that the information Cell (method ReadFieldCell) is repeated excessively (indeed the field "CellField" was created on cells), but remember that in the MED file standard, nothing prevents a field from being based on cells but also simultaneously on nodes, or Gauss points ...</div>
#
# Compare `cellField_read` and `cellField0`. Problem, due to the constraint on MED file numbering, we have lost the original numbering. Or more precisely, there is no standard way to retrieve the original numbering. So a `MEDCouplingFieldDouble.isEqual()` is not enough. Let's use a `MEDCouplingFieldDouble.substractInPlaceDM()` which performs for us a renumbering following a particular policy (policy, see html doc). To do this, make a deep copy of `cellField` to `cellFieldCpy` and perform a `substractInPlaceDM` (DM for "Different Meshes") on this copy (unlike `substract` which only works if they share the same mesh):

cellFieldCpy = cellField.deepCopy()
cellFieldCpy.substractInPlaceDM(cellField_read, 10, 1e-12)
cellFieldCpy.getArray().abs()
print(cellFieldCpy.getArray().isUniform(0.0, 1e-12))

# Perform the same work on "NodeField" as done earlier on "CellField". The difference here is that there will be duplication of information at the boundary, because the boundary nodes are shared on both sides:

nodeField0_read = mc.ReadFieldNode("proc0.med", "mesh", 0, "NodeField", 5, 6)
nodeField1_read = mc.ReadFieldNode("proc1.med", "mesh", 0, "NodeField", 5, 6)
nodeField_read = mc.MEDCouplingFieldDouble.MergeFields(
    [nodeField0_read, nodeField1_read]
)

# <div class="alert alert-block alert-success">
# <b>Note:</b> In this part, we have re-read the mesh a second time, which can be penalizing ...</div>
#
# Invoke `MEDCouplingUMesh.mergeNodes()` on `nodeField_read` to remove the duplicated nodes. Make a deep copy called `nodeFieldCpy` of `nodeField` and again remove duplicates:

nodeField_read.mergeNodes(1e-10)
nodeFieldCpy = nodeField.deepCopy()
nodeFieldCpy.mergeNodes(1e-10)

# <div class="alert alert-block alert-success">
# <b>Note:</b> Note that mergeNodes() has two precision parameters (epsilons), the first, classic, on the absolute distance between nodes, and the other on the tolerance accepted on the values of the field. If the value of the field of two nodes to be merged exceeds this second epsilon, an exception is raised.</div>
#
# Compare `nodeFieldCpy` and `nodeField_read` again using `MEDCouplingFieldDouble.substractInPlaceDM()`:

nodeFieldCpy.substractInPlaceDM(nodeField_read, 10, 1e-12)
print(nodeFieldCpy.getArray().isUniform(0.0, 1e-12))


# ## Reading and Merging the 2 Separate MED Files (less easy, but more optimal)
#
# Here, we need to perform a more systematic and potentially more general method of file merging. For large files, this approach is preferable. Besides performance, this approach has the advantage of being able to add information.
#
# With the advanced API, read the meshes of the two files "proc0.med" and "proc1.med" and aggregate the result in an instance `mergeMLMesh` of `MEDFileUMesh`. Process all dimension levels (even if here there is only one) using the method `MEDFileUMesh.getNonEmptyLevels()` on the instance coming from "proc0.med".
#
# The solution given below is as generic as possible, as it also handles different timesteps and different geometric types:


# +
def load_mesh_and_fields(fileNames):
    msML = [mc.MEDFileMesh.New(fname) for fname in fileNames]
    fsML = [mc.MEDFileFields.New(fname) for fname in fileNames]
    return msML, fsML


def merge_meshes(msML):
    mergeMLMesh = mc.MEDFileUMesh()
    o2nML = {}  # Initialize o2nML as a dictionary to store mappings for each level
    for lev in msML[0].getNonEmptyLevels():
        cs = [mML.getCoords() for mML in msML]
        mergeMLMesh.setCoords(mc.DataArrayDouble.Aggregate(cs))
        ms = [mML.getMeshAtLevel(lev) for mML in msML]
        m = mc.MEDCouplingUMesh.MergeUMeshes(ms)
        m.setCoords(mergeMLMesh.getCoords())
        o2nML[lev] = m.sortCellsInMEDFileFrmt()  # Store mapping for the current level
        mergeMLMesh.setMeshAtLevel(lev, m)
    return mergeMLMesh, o2nML


def merge_fields(fsML, mergeMLMesh, o2nML):
    mergeMLFields = mc.MEDFileFields()
    for fieldName in fsML[0].getFieldsNames():
        fmts = [fML[fieldName] for fML in fsML]
        mergeField = mc.MEDFileFieldMultiTS()
        for dt, it, tim in fmts[0].getTimeSteps():
            fts = [fmt[dt, it] for fmt in fmts]
            for typp in fts[0].getTypesOfFieldAvailable():
                arr1s = aggregate_data(fts, typp)
                for lev in o2nML:  # Ensure we have mappings for all levels
                    arr = mc.DataArrayDouble.Aggregate(arr1s)
                    if typp == mc.ON_CELLS and lev in o2nML:
                        arr.renumberInPlace(
                            o2nML[lev]
                        )  # Use the correct mapping for each level
                    mcf = mc.MEDCouplingFieldDouble(typp, mc.ONE_TIME)
                    mcf.setName(fieldName)
                    mcf.setTime(tim, dt, it)
                    mcf.setArray(arr)
                    mcf.setMesh(mergeMLMesh.getMeshAtLevel(lev))
                    mcf.checkConsistencyLight()
                    mergeField.appendFieldNoProfileSBT(mcf)
        mergeMLFields.pushField(mergeField)
    return mergeMLFields


def aggregate_data(fts, typp):
    arr1s = []
    for ft in fts:
        for geoTyp, smth in ft.getFieldSplitedByType():
            if geoTyp != mc.NORM_ERROR:
                smth1 = [elt for elt in smth if elt[0] == mc.ON_CELLS]
                arr2s = [
                    ft.getUndergroundDataArray()[elt[1][0] : elt[1][1]] for elt in smth1
                ]
                arr1s.append(mc.DataArrayDouble.Aggregate(arr2s))
            else:
                smth = [
                    elt for elt in ft.getFieldSplitedByType() if elt[0] == mc.NORM_ERROR
                ]
                arr2 = mc.DataArrayDouble.Aggregate(
                    [
                        ft.getUndergroundDataArray()[elt[1][0][1][0] : elt[1][0][1][1]]
                        for elt in smth
                    ]
                )
                arr1s.append(arr2)
    return arr1s


def write_merged_data(mergeMLMesh, mergeMLFields, fileName="merge.med"):
    mergeMLMesh.write(fileName, 2)
    mergeMLFields.write(fileName, 0)


# Main execution
fileNames = ["proc0.med", "proc1.med"]
msML, fsML = load_mesh_and_fields(fileNames)
mergeMLMesh, o2nML = merge_meshes(msML)
mergeMLFields = merge_fields(fsML, mergeMLMesh, o2nML)
write_merged_data(mergeMLMesh, mergeMLFields)
