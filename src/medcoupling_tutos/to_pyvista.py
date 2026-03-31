import medcoupling as mc
import numpy as np
import pyvista as pv

MC_TO_PV_CELLTYPE = {
    mc.NORM_POINT1: pv.CellType.VERTEX,
    mc.NORM_QUAD4: pv.CellType.QUAD,
    mc.NORM_TRI3: pv.CellType.TRIANGLE,
    mc.NORM_HEXA8: pv.CellType.HEXAHEDRON,
    mc.NORM_SEG2: pv.CellType.LINE,
    mc.NORM_SEG3: pv.CellType.QUADRATIC_EDGE,
    mc.NORM_TETRA4: pv.CellType.TETRA,
    mc.NORM_POLYGON: pv.CellType.POLYGON,
    mc.NORM_QPOLYG: pv.CellType.QUADRATIC_POLYGON,
    mc.NORM_POLYHED: pv.CellType.POLYHEDRON,
}

MC_DIM = {
    mc.NORM_POINT1: 0,
    mc.NORM_QUAD4: 2,
    mc.NORM_TRI3: 2,
    mc.NORM_HEXA8: 3,
    mc.NORM_SEG2: 1,
    mc.NORM_SEG3: 1,
    mc.NORM_TETRA4: 3,
    mc.NORM_POLYGON: 2,
    mc.NORM_QPOLYG: 2,
    mc.NORM_POLYHED: 3,
}


def _mc_ph_to_vtk(cell_co: np.ndarray) -> np.ndarray:
    seps = np.flatnonzero(cell_co == -1)
    seps = np.r_[0, seps, len(cell_co)]
    faces = [
        np.r_[f_end - f_start - 1, cell_co[f_start + 1 : f_end]]
        for f_start, f_end in zip(seps[:-1], seps[1:])
    ]

    vtk_co = np.r_[*faces]
    assert sum(f[0] for f in faces) + len(faces) == len(vtk_co)

    res = np.r_[len(vtk_co) + 1, len(faces), vtk_co]
    return res


def _to_unstructured(
    mesh: mc.MEDCouplingUMesh, coords: np.ndarray
) -> pv.UnstructuredGrid:
    offsets = mesh.getNodalConnectivityIndex().toNumPyArray()
    cell_length = offsets[1:] - offsets[:-1] - 1

    cell_types_idx = offsets[:-1]
    connectivity = np.array(mesh.getNodalConnectivity().toNumPyArray())
    mc_cell_types = np.array(connectivity[cell_types_idx])

    # Split off polyhedrons
    # NOTE: polyhedrons should be higher types id
    # TODO: I should check that and isolate the POLYHED case...
    ind_l = np.searchsorted(mc_cell_types, mc.NORM_POLYHED, side="left")
    ind_r = np.searchsorted(mc_cell_types, mc.NORM_POLYHED, side="right")
    cell_types = np.array(mc_cell_types)

    # Non polyhedron case
    types_idx = dict()
    for mc_type in MC_TO_PV_CELLTYPE:
        types_idx[mc_type] = cell_types == mc_type
    for mc_type, pv_type in MC_TO_PV_CELLTYPE.items():
        cell_types[types_idx[mc_type]] = pv_type

    connectivity_npl = connectivity[: offsets[ind_l]]
    connectivity_npl[cell_types_idx[:ind_l]] = cell_length[:ind_l]
    connectivity_npr = connectivity[offsets[ind_r] :]
    connectivity_npr[cell_types_idx[ind_r:]] = cell_length[ind_r:]

    # Polyhedron case
    phs = [
        _mc_ph_to_vtk(connectivity[cell_start:cell_end])
        for cell_start, cell_end in zip(
            offsets[ind_l:ind_r], offsets[ind_l + 1 : ind_r + 1]
        )
    ]
    if len(phs) > 0:
        connectivity_ph = np.r_[*phs]
    else:
        connectivity_ph = np.array([], dtype=int)

    # Combination
    connectivity = np.r_[connectivity_npl, connectivity_npr, connectivity_ph]

    pv_mesh = pv.UnstructuredGrid(connectivity, cell_types, coords)
    return pv_mesh


def _to_polydata(mesh: mc.MEDCouplingUMesh, coords: np.ndarray) -> pv.PolyData:
    offsets = mesh.getNodalConnectivityIndex().toNumPyArray()
    cell_length = offsets[1:] - offsets[:-1] - 1
    cell_types_idx = offsets[:-1]

    connectivity = np.array(mesh.getNodalConnectivity().toNumPyArray())
    any_type = connectivity[0]
    connectivity[cell_types_idx] = cell_length

    if MC_DIM[any_type] == 0:
        return pv.PolyData(coords, verts=connectivity)
    if MC_DIM[any_type] == 1:
        return pv.PolyData(coords, lines=connectivity)
    if MC_DIM[any_type] == 2:
        return pv.PolyData(coords, faces=connectivity)
    else:
        raise Exception(f"{any_type} if not in known types: {MC_DIM=}")


def to_pv(
    mesh: mc.MEDCouplingUMesh | mc.MEDCouplingFieldDouble,
) -> pv.UnstructuredGrid | pv.PolyData:
    if isinstance(mesh, mc.MEDCouplingFieldDouble):
        field: mc.MEDCouplingFieldDouble = mesh
        mesh = mesh.getMesh()
    else:
        field = None
        fname = None

    # Prepare mesh : make copy and sort cell types
    mesh = mesh.deepCopyConnectivityOnly()
    permut = mesh.sortCellsInMEDFileFrmt()

    coords = np.array(mesh.getCoords().toNumPyArray())
    if coords.shape[1] == 3:
        pv_mesh = _to_unstructured(mesh, coords)
    elif coords.shape[1] < 3:
        coords = np.c_[coords, np.zeros((coords.shape[0], 3 - coords.shape[1]))]
        pv_mesh = _to_polydata(mesh, coords)
    else:
        raise Exception(f"The coords shape is not valid: {coords.shape=}")

    if field is not None:
        fname = mesh.getName()
        farr = field.getArray().toNumPyArray()
        if fname is None or fname == "":
            fname = "Field"
        if field.getTypeOfField() == mc.ON_CELLS:
            permut = permut.toNumPyArray()
            p = np.empty_like(permut)
            p[permut] = np.arange(permut.size)
            pv_mesh.cell_data[fname] = farr[p]
        elif field.getTypeOfField() == mc.ON_NODES:
            pv_mesh.point_data[fname] = farr
        else:
            raise NotImplementedError(f"{field.getTypeOfField()=} is not supported.")

    return pv_mesh
