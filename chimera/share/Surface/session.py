# -----------------------------------------------------------------------------
# Save surface pieces having attribute save_in_session = True.
#
def save_surface_state(file):

  slist = []
  from chimera import openModels as om
  from _surface import SurfaceModel
  for s in om.list(modelTypes = [SurfaceModel]):
    plist = [p for p in s.surfacePieces
             if hasattr(p, 'save_in_session') and p.save_in_session]
    if plist:
      slist.append((s,plist)) 
  if len(slist) == 0:
    return

  smslist = []
  for s,plist in slist:
    sms = Surface_Model_State()
    sms.state_from_surface_model(s,plist)
    smslist.append(sms)
  
  from SessionUtil import objecttree
  t = objecttree.instance_tree_to_basic_tree(smslist)

  file.write('\n')
  file.write('def restore_surfaces():\n')
  file.write(' surface_state = \\\n')
  objecttree.write_basic_tree(t, file, indent = '  ')
  file.write('\n')
  file.write(' from Surface import session\n')
  file.write(' session.restore_surface_state(surface_state)\n')
  file.write('\n')
  file.write('try:\n')
  file.write(' restore_surfaces()\n')
  file.write('except:\n')
  file.write(" reportRestoreError('Error restoring surfaces')\n")
  file.write('\n')
  
# -----------------------------------------------------------------------------
#
def restore_surface_state(surface_state):

  from SessionUtil.stateclasses import Model_State, Xform_State, Arrays_State

  classes = (
    Surface_Model_State,
    Surface_Piece_State,
    Model_State,
    Xform_State,
    Arrays_State,
    )
  name_to_class = {}
  for c in classes:
    name_to_class[c.__name__] = c

  from SessionUtil import objecttree
  stlist = objecttree.basic_tree_to_instance_tree(surface_state, name_to_class)

  for st in stlist:
    s = st.create_object()
    for p in s.surfacePieces:
      p.save_in_session = True

# -----------------------------------------------------------------------------
#
class Surface_Model_State:
  
  version = 3

  state_attributes = ('name',
		      'surface_model',
                      'arrays',
                      'pieces_are_selectable',
                      'surface_piece_states',
		      'version',
		      )

  # ---------------------------------------------------------------------------
  #
  def state_from_surface_model(self, sm, plist):

    self.name = sm.name

    from SessionUtil.stateclasses import Model_State, Arrays_State
    ms = Model_State()
    ms.state_from_model(sm)
    self.surface_model = ms
    self.arrays = Arrays_State()
    
    self.pieces_are_selectable = sm.piecesAreSelectable

    splist = []
    for p in plist:
      sp = Surface_Piece_State()
      sp.state_from_surface_piece(p, self.arrays)
      splist.append(sp)
    self.surface_piece_states = splist

  # ---------------------------------------------------------------------------
  #
  def create_object(self):

    from _surface import SurfaceModel
    sm = SurfaceModel()
    arrays = self.arrays if self.version >= 3 else None
    
    sm.piecesAreSelectable = self.pieces_are_selectable

    for ps in self.surface_piece_states:
      ps.create_object(sm, arrays)

    sms = self.surface_model
    from SimpleSession import modelOffset
    from chimera import openModels as om
    om.add([sm], baseId = sms.id + modelOffset, subid = sms.subid)
    sms.restore_state(sm)

    return sm

# -----------------------------------------------------------------------------
#
class Surface_Piece_State:
  
  version = 4

  state_attributes = ('vertices',
		      'triangles',
                      'normals',
                      'placement',
                      'color',
                      'per_vertex_colors',
                      'triangle_and_edge_mask',
                      'displayStyle',
                      'useLighting',
                      'twoSidedLighting',
                      'lineThickness',
                      'dotSize',
                      'smoothLines',
                      'transparencyBlendMode',
		      'version',
		      )

  # ---------------------------------------------------------------------------
  #
  def state_from_surface_piece(self, p, arrays):

    if p.placement is None:
      varray, tarray = p.geometry
      narray = p.normals
    else:
      varray, narray = p.unplacedVerticesAndNormals
      tarray = p.geometry[1]
    
    from SessionUtil import float32, int32
    self.vertices = arrays.id(varray, float32)
    self.triangles = arrays.id(tarray, int32)
    self.normals = arrays.id(narray, float32)
    self.placement = p.placement
    self.color = p.color
    self.per_vertex_colors = arrays.id(p.vertexColors, float32)
    self.triangle_and_edge_mask = arrays.id(p.triangleAndEdgeMask, int32)
    for attr in ('displayStyle', 'useLighting', 'twoSidedLighting',
                 'lineThickness', 'dotSize', 'smoothLines',
                 'transparencyBlendMode'):
      setattr(self, attr, getattr(p, attr))

  # ---------------------------------------------------------------------------
  #
  def create_object(self, surf_model, arrays):

    if self.version >= 4:
      varray = arrays.array(self.vertices)
      tarray = arrays.array(self.triangles)
    elif self.version >= 3:
      from SessionUtil import string_to_array, float32, int32
      varray = string_to_array(self.vertices, float32, 3)
      tarray = string_to_array(self.triangles, int32, 3)
    else:
      varray = self.vertices
      tarray = self.triangles
    p = surf_model.addPiece(varray, tarray, self.color)

    if self.version >= 4 and not self.placement is None:
      p.placement = self.placement

    if self.version >= 4:
      rgba = arrays.array(self.per_vertex_colors)
    elif self.version >= 3:
      rgba = string_to_array(self.per_vertex_colors, float32, 4)
    else:
      rgba = self.per_vertex_colors
    p.vertexColors = rgba

    if self.version >= 4:
      narray = arrays.array(self.normals)
      # Avoid error caused by old sessions where normals not saved correctly in session.
      if len(narray) == len(varray):
        p.normals = narray
      p.triangleAndEdgeMask = arrays.array(self.triangle_and_edge_mask)
    elif self.version >= 3:
      p.normals = string_to_array(self.normals, float32, 3)
      p.triangleAndEdgeMask = string_to_array(self.triangle_and_edge_mask,int32)
    elif self.version >= 2:
      p.normals = self.normals
      p.triangleAndEdgeMask = self.triangle_and_edge_mask
      
    if self.version >= 2:
      for attr in ('displayStyle', 'useLighting', 'twoSidedLighting',
                   'lineThickness', 'dotSize', 'smoothLines',
                   'transparencyBlendMode'):
        setattr(p, attr, getattr(self, attr))

    return p
