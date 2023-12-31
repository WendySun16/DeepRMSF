# -----------------------------------------------------------------------------
# Produce hexagonal lattices on an icosahedron.  The hexagons are bent where
# they cross the edges of the icosahedron.
#
# These lattices are described at
#
#       http://viperdb.scripps.edu/icos_server.php?icspage=paradigm
#
def show_hk_lattice(h, k, radius, orientation = '222',
                    color_rgba = (1,1,1,1), mesh_line_thickness = 1,
                    sphere_factor = 0, replace = True):

    varray, tarray, hex_edges = hk_icosahedron_lattice(h,k,radius,orientation)
    interpolate_with_sphere(varray, radius, sphere_factor)

    # Make new surface model or find an existing one.
    sm = None
    open = False
    from chimera import openModels
    from _surface import SurfaceModel
    if replace:
        mlist = openModels.list(modelTypes = [SurfaceModel])
        mlist = [m for m in mlist if hasattr(m, 'hkcage')]
        if mlist:
            sm = mlist[0]
            for p in sm.surfacePieces:
                if hasattr(p, 'hkcage'):
                    sm.removePiece(p)
            open = True
    if sm is None:
        sm = SurfaceModel()
        sm.name = 'h = %d, k = %d lattice' % (h,k)
        sm.hkcage = True

    p = sm.addPiece(varray, tarray, color_rgba)
    p.hkcage = True
    p.lineThickness = mesh_line_thickness
    p.displayStyle = p.Mesh
    p.setEdgeMask(hex_edges)    # Hide spokes of hexagons.

    if not open:
        openModels.add([sm])

    return sm

# -----------------------------------------------------------------------------
#
def hk_icosahedron_lattice(h, k, radius, orientation):

    # Find triangles for the hk lattice covering one asymmetric unit equilateral triangle.
    # The asym unit triangle (corners) and hk lattice triangles are in the xy plane in 3-d.
    corners, triangles, t_hex_edges = hk_triangle(h, k)

    from Icosahedron import icosahedron_geometry
    ivarray, itarray = icosahedron_geometry(orientation)

    # Map the 2d hk asymmetric unit triangles onto each face of an icosahedron
    tlist = []
    for i0,i1,i2 in itarray:
        face = ivarray[i0], ivarray[i1], ivarray[i2]
        tmap = triangle_map(corners, face)
        tlist.extend(map_triangles(tmap, triangles))

    # Convert from triangles defined by 3 vertex points, to an array of
    # unique vertices and triangles as 3 indices into the unique vertex list.
    va, ta = surface_geometry(tlist, tolerance = 1e-5)

    # Scale to requested radius
    from numpy import multiply, array, intc
    multiply(va, radius, va)
    
    # Compute the edge mask to show just the hexagon edges.
    hex_edges = array(t_hex_edges * len(itarray), intc)
    
    return va, ta, hex_edges

# -----------------------------------------------------------------------------
# Calculate a list of 2-d triangles that is the part of the hk lattice within
# a single asymmetric unit equilateral triangle in the plane.  Triangles that
# would cross the boundary of the asymmetric unit triangle are clipped.
#
# The returned corners are the 3 vertices of the asymmetric unit triangle.
# The returned asym unit and hk triangles are in the xy plane in 3 dimensions (z=0).
# A returned edge mask indicates which of the hk lattice triangle edges
# are hexagon edges.
#
def hk_triangle(h, k):

    # Multiply h,k by 3 so hexagon corners have integer coordinates for
    # exact intersection calculations.
    corners2d = ((0,0), (3*h,3*k), (-3*k,3*(h+k)))

    hex_corner_offset = ((2,-1), (1,1), (-1,2), (-2,1), (-1,-1), (1,-2))
    kmax = max(k, h+k)
    triangles2d = []
    hex_edges = []
    for k0 in range(kmax+1):
        for h0 in range(-k0, h+1):
            for c in range(6):
                h1o, k1o = hex_corner_offset[c]
                h2o, k2o = hex_corner_offset[(c+1)%6]
                tri = ((3*h0,3*k0), (3*h0+h1o,3*k0+k1o), (3*h0+h2o,3*k0+k2o))
                ti, he = triangle_intersection(tri, corners2d, 2)
                triangles2d.extend(ti)
                hex_edges.extend(he)

    corners = hk3_to_xyz(corners2d)
    triangles = map(lambda t: hk3_to_xyz(t), triangles2d)
    
    return corners, triangles, hex_edges
    
# -----------------------------------------------------------------------------
# Triangulate the portion of triangle t1 inside t2.  The triangles are specified
# by 3 vertex points and are in 2 dimensions.  Only the cases that occur in the
# hk icosahedral grids are handled.
#
def triangle_intersection(t1, t2, edge_mask):

    interior_vertices = []
    exterior_vertices = []
    boundary_vertices = []
    locations = {-1:exterior_vertices, 0:boundary_vertices, 1:interior_vertices}
    for k in range(3):
        loc = vertex_in_triangle(t1[k], t2)
        locations[loc].append(k)

    iv = len(interior_vertices)
    ev = len(exterior_vertices)
    if iv == 0 and ev > 0:
        return [], []
    if ev == 0:
        return [t1], [edge_mask]

    # Have at least one exterior and one interior vertex.  Need new triangles.
    if iv == 1 and ev == 1:
        i = interior_vertices[0]
        b = boundary_vertices[0]
        e = exterior_vertices[0]
        thalf = list(t1)
        thalf[e] = cut_point(t1[i], t1[e], t2)
        return [thalf], [mask_edge(edge_mask, (e,b))]
    
    if iv == 1 and ev == 2:
        i = interior_vertices[0]
        e1 = exterior_vertices[0]
        e2 = exterior_vertices[1]
        tpiece = list(t1)
        tpiece[e1] = cut_point(t1[i], t1[e1], t2)
        tpiece[e2] = cut_point(t1[i], t1[e2], t2)
        return [tpiece], [mask_edge(edge_mask, (e1,e2))]

    if iv == 2 and ev == 1:
        i1 = interior_vertices[0]
        i2 = interior_vertices[1]
        e = exterior_vertices[0]
        tpiece1 = list(t1)
        tpiece1[e] = cut_point(t1[i1], t1[e], t2)
        em1 = mask_edge(edge_mask, (i2,e))
        tpiece2 = list(t1)
        tpiece2[e] = cut_point(t1[i2], t1[e], t2)
        tpiece2[i1] = tpiece1[e]
        em2 = mask_edge(edge_mask, (i1,e), (i1,i2))
        return [tpiece1, tpiece2], [em1, em2]
    

    raise ValueError, 'Unexpected triangle intersection'

# -----------------------------------------------------------------------------
# Vertex and triangle are in two dimensions with triangle defined by 3 corner
# vertices.
#
def vertex_in_triangle(v, t):

    v0, v1, v2 = t
    v01p = (-v1[1]+v0[1], v1[0]-v0[0])
    v12p = (-v2[1]+v1[1], v2[0]-v1[0])
    v20p = (-v0[1]+v2[1], v0[0]-v2[0])
    from numpy import subtract, dot as inner_product
    sides = (inner_product(subtract(v,v0), v01p),
             inner_product(subtract(v,v1), v12p),
             inner_product(subtract(v,v2), v20p))
    inside = len(filter(lambda s: s > 0, sides))
    outside = len(filter(lambda s: s < 0, sides))
    if inside == 3:
        return 1     # Interior
    elif outside == 0:
        return 0     # Boundary
    
    return -1        # Outside

# -----------------------------------------------------------------------------
# Find intersection of segment u->v with triangle boundary.
# u and v must be an interior and an exterior point.
#
def cut_point(u, v, tri):

    for e in range(3):
        ip = segment_intersection(u, v, tri[e], tri[(e+1)%3])
        if ip:
            return ip
    raise ValueError, 'No intersection %s %s %s' % (u, v, tri)
    return None

# -----------------------------------------------------------------------------
# Find intersection of segment ab with segment cd.
# Return (x,y) or None if no intersectin.
#
def segment_intersection(a, b, c, d):

    m11 = b[0]-a[0]
    m21 = b[1]-a[1]
    m12 = c[0]-d[0]
    m22 = c[1]-d[1]
    det = m11*m22 - m12*m21
    if det == 0:
        return None

    y1 = c[0]-a[0]
    y2 = c[1]-a[1]
    f = float(m22*y1 - m12*y2)/det
    g = float(-m21*y1 + m11*y2)/det
    if f < 0 or f > 1 or g < 0 or g > 1:
        return None

    p = (a[0] + m11*f, a[1] + m21*f)
    return p

# -----------------------------------------------------------------------------
# Mask out edges given by pair of vertex indices (0-2).  Bits 0, 1, and 2
# correspond to edges 0-1, 1-2, and 2-0 respectively.
#
def mask_edge(edge_mask, *edges):

    ebits = {(0,1):1, (1,0):1, (1,2):2, (2,1):2, (2,0):4, (0,2):4}
    emask = edge_mask
    for e in edges:
        emask &= ~ebits[e]
    return emask

# -----------------------------------------------------------------------------
# Shear transform 2d hk points to points on the xy plane in 3 dimensions (z=0).
#
def hk3_to_xyz(hklist):

    from math import sqrt
    hx = sqrt(3)/6
    hy = 0.5/3
    ky = 1.0/3
    xyz_list = map(lambda hk: (hk[0]*hx,hk[1]*ky+hk[0]*hy,0), hklist)
    return xyz_list

# -----------------------------------------------------------------------------
# Compute the 3 by 4 transform matrix mapping one 3-d triangle to another.
#
def triangle_map(tri1, tri2):

    from numpy import zeros, subtract, dot as matrix_multiply, float

    f1 = zeros((3,3), float)
    f1[:,0], f1[:,1] = subtract(tri1[1], tri1[0]), subtract(tri1[2], tri1[0])
    f1[:,2] = cross_product(f1[:,0], f1[:,1])

    f2 = zeros((3,3), float)
    f2[:,0], f2[:,1] = subtract(tri2[1], tri2[0]), subtract(tri2[2], tri2[0])
    f2[:,2] = cross_product(f2[:,0], f2[:,1])
    
    from numpy.linalg import inv as inverse
    f1inv = inverse(f1)

    tmap = zeros((3,4), float)
    tmap[:,:3] = matrix_multiply(f2, f1inv)
    tmap[:,3] = subtract(tri2[0], matrix_multiply(tmap[:,:3],tri1[0]))

    return tmap
    
# -----------------------------------------------------------------------------
# Apply a 3x4 affine transformation to vertices of triangles.
#
def map_triangles(tmap, triangles):

    from Matrix import apply_matrix
    tri = map(lambda t: map(lambda v: apply_matrix(tmap, v), t), triangles)
    return tri
    
# -----------------------------------------------------------------------------
#
def cross_product(u, v):

    return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
    
# -----------------------------------------------------------------------------
# Take a list of triangles where each triangle is specified by 3 xyz vertex
# positions and convert to a vertex and triangle array where the triangle
# array contains indices into the vertex array.  Vertices in the original
# triangle data that are close (within tolerance) are merged into a single
# vertex.
#
def surface_geometry(triangles, tolerance = 1e-5):

    from numpy import array, reshape, single as floatc, intc
    varray = reshape(triangles, (3*len(triangles),3)).astype(floatc)
    
    uindex = {}
    unique = []
    from _closepoints import find_close_points, BOXES_METHOD
    for v in range(len(varray)):
        if not v in uindex:
            i1, i2 = find_close_points(BOXES_METHOD, varray[v:v+1,:], varray,
                                       tolerance)
            for i in i2:
                if not i in uindex:
                    uindex[i] = len(unique)
            unique.append(varray[v])

    uvarray = array(unique, floatc)
    tlist = map(lambda t: (uindex[3*t],uindex[3*t+1],uindex[3*t+2]),
                range(len(triangles)))
    tarray = array(tlist, intc)
    
    return uvarray, tarray

# -----------------------------------------------------------------------------
# Radially interpolate vertex points a certain factor towards a sphere of
# given radius.
#
def interpolate_with_sphere(varray, radius, sphere_factor):

    if sphere_factor == 0:
        return

    from math import sqrt
    for v in range(len(varray)):
        x,y,z = varray[v]
        r = sqrt(x*x + y*y + z*z)
        if r > 0:
            ri = r * (1-sphere_factor) + radius*sphere_factor
            f = ri/r
            varray[v,:] = (f*x, f*y, f*z)
