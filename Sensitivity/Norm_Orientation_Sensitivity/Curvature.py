'''
Investigating the spatial interaction of immune cells in colon cancer

Curvature: Calculates Curvature and Normal vectors of a given domain.
            Courtesy of http://jsdokken.com/

Modified by: Navid Mohammad Mirzaei https://sites.google.com/view/nmirzaei
                               https://github.com/nmirzaei

(c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
'''
###############################################################
#Importing required functions
###############################################################
import numpy
# import dolfin
from dolfin import *
from dolfin_adjoint import *
###############################################################
def Curvature(mesh):

    ###############################################################
    #Returns a the boundary representation of the CG-1 function v
    ###############################################################
    def mesh_to_boundary(v, b_mesh):
        ###############################################################
        # Extract the underlying volume and boundary meshes
        ###############################################################
        mesh = v.function_space().mesh()
        ###############################################################

        ###############################################################
        # We use a Dof->Vertex mapping to create a global
        # array with all DOF values ordered by mesh vertices
        ###############################################################
        DofToVert = dof_to_vertex_map(v.function_space())
        VGlobal = numpy.zeros(v.vector().size())
        ###############################################################

        ###############################################################
        #Creating a map compatible with parallel computation
        ###############################################################
        vec = v.vector().get_local()
        for i in range(len(vec)):
            Vert = MeshEntity(mesh, 0, DofToVert[i])
            globalIndex = Vert.global_index()
            VGlobal[globalIndex] = vec[i]
        VGlobal = SyncSum(VGlobal)
        ###############################################################

        ###############################################################
        # Use the inverse mapping to se the DOF values of a boundary
        # function
        ###############################################################
        surface_space = FunctionSpace(b_mesh, "CG", 1)
        surface_function = Function(surface_space)
        mapa = b_mesh.entity_map(0)
        DofToVert = dof_to_vertex_map(FunctionSpace(b_mesh, "CG", 1))

        LocValues = surface_function.vector().get_local()
        for i in range(len(LocValues)):
            VolVert = MeshEntity(mesh, 0, mapa[int(DofToVert[i])])
            GlobalIndex = VolVert.global_index()
            LocValues[i] = VGlobal[GlobalIndex]

        surface_function.vector().set_local(LocValues)
        surface_function.vector().apply('')
        return surface_function
        ###############################################################

    ###############################################################
    #Map bulk values to boundary
    ###############################################################
    def vector_mesh_to_boundary(func, b_mesh):
        v_split = func.split(deepcopy=True)
        v_b = []
        for v in v_split:
            v_b.append(mesh_to_boundary(v, b_mesh))
        Vb = VectorFunctionSpace(b_mesh, "CG", 1)
        vb_out = Function(Vb)
        scalar_to_vec = FunctionAssigner(Vb, [v.function_space() for
                                                      v in v_b])
        scalar_to_vec.assign(vb_out, v_b)
        return vb_out
    ###############################################################

    ###############################################################
    #Returns sum of vec over all mpi processes.
    #Each vec vector must have the same dimension for each MPI process
    ###############################################################
    def SyncSum(vec):
        comm = MPI.comm_world
        NormalsAllProcs = numpy.zeros(comm.Get_size() * len(vec), dtype=vec.dtype)
        comm.Allgather(vec, NormalsAllProcs)

        out = numpy.zeros(len(vec))
        for j in range(comm.Get_size()):
            out += NormalsAllProcs[len(vec) * j:len(vec) * (j + 1)]
        return out
    ###############################################################

    ###############################################################
    #Mapping boundary values to mesh values
    ###############################################################
    def boundary_to_mesh(f, mesh):
        b_mesh = f.function_space().mesh()
        SpaceV = FunctionSpace(mesh, "CG", 1)
        SpaceB = FunctionSpace(b_mesh, "CG", 1)

        F = Function(SpaceV)
        GValues = numpy.zeros(F.vector().size())

        map = b_mesh.entity_map(0)  # Vertex map from boundary mesh to parent mesh
        d2v = dof_to_vertex_map(SpaceB)
        v2d = vertex_to_dof_map(SpaceV)

        dof = SpaceV.dofmap()
        imin, imax = dof.ownership_range()

        for i in range(f.vector().local_size()):
            GVertID = Vertex(b_mesh, d2v[i]).index()  # Local Vertex ID for given dof on boundary mesh
            PVertID = map[GVertID]  # Local Vertex ID of parent mesh
            PDof = v2d[PVertID]  # Dof on parent mesh
            value = f.vector()[i]  # Value on local processor
            GValues[dof.local_to_global_index(PDof)] = value
        GValues = SyncSum(GValues)

        F.vector().set_local(GValues[imin:imax])
        F.vector().apply("")
        return F
     ###############################################################


    ###############################################################
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)
    ###############################################################

    ###############################################################
    A.ident_zeros()
    nh = Function(V)
    ###############################################################

    ###############################################################
    solve(A, nh.vector(), L)
    File("Results/nh.pvd") << nh
    ###############################################################

    ###############################################################
    bmesh = BoundaryMesh(mesh, "exterior")
    nb = vector_mesh_to_boundary(nh, bmesh)
    Q = FunctionSpace(bmesh, "P", 1)
    ###############################################################

    ###############################################################
    p, q = TrialFunction(Q), TestFunction(Q)
    a = inner(p,q)*dx
    l = inner(0.5*div(nb), q)*dx
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)
    A.ident_zeros()
    kappab = Function(Q)
    solve(A, kappab.vector(), L)
    kappa = boundary_to_mesh(kappab, mesh)
    File("test/kappa.pvd") << kappa
    return kappa, nh
    ###############################################################
###############################################################
