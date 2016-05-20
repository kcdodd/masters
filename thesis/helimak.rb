

# Define all subroutines first
# Make arc in xy plane with perpendicular ends (for extruding)
def normal_arc(thBeg, thEnd, radius, numseg)
  # thBeg, thEnd  =  begin & end angles with x-axis (degrees)
  
  # Calc various angles in radians
  thTot  =  thEnd-thBeg                # Total arc angle
  dTh    =  thTot.degrees/(numseg-1)   # full angle increment per segment
  
  radius_reduced = radius*Math.cos(dTh/2.0) # for first/last point must have slightly smaller radius

  pts  =  []
  theta  = thBeg.degrees
  
  pts << Geom::Point3d.new(radius_reduced*Math.cos(theta), radius_reduced*Math.sin(theta),0)
  
  theta = theta + dTh/2.0

  for i in 1..(numseg-1) do
	
	pts <<  Geom::Point3d.new(radius*Math.cos(theta), radius*Math.sin(theta),0)
	theta = theta + dTh
  end
  
  theta = theta - dTh/2.0
  
  pts <<  Geom::Point3d.new(radius_reduced*Math.cos(theta), radius_reduced*Math.sin(theta),0)

  return pts
end

# takes points and extrudes it in an arc about the z axis, begining at arbitrary thBeg
# assumes face starts in x-z plane
def extrude_arc(ent, pts, thBeg, thEnd, soften_angle)
	gr = ent.add_group
	
	grent = gr.entities

	thTot = thEnd-thBeg
	
	curve = grent.add_curve pts
	
	face = grent.add_face curve

	# extrude starting at x axis
	extr = grent.add_curve normal_arc(0, thTot, 100, 51)
	
	# Extrude plasma face into torus
	arc = face.followme extr
	
	grent.erase_entities extr
	
	#now rotate whole extrusion about the z-axis to the correct starting position
	trans = Geom::Transformation.rotation [0, 0, 0], [0, 0, 1], thBeg.degrees
	
	gr.transformation = trans
	
	soften(gr, soften_angle, soften_angle)
	
	return gr
end

# performs soft/smooth function for a group and all sub-groups
def soften(group, angleMin, angleMax)

	group.entities.each{|e|

		if (e.class == Sketchup::Group)
			# appply recursivly
			soften(e, angleMin, angleMax)
			
		elsif (e.class == Sketchup::Edge)
		
			if (e.faces.length == 2)
				# only test if number of faces is at least 2
				test_angle = e.faces[0].normal.angle_between e.faces[1].normal

				if (test_angle <= angleMin.degrees)
				
					e.soft = true
					e.smooth = true
					
				elsif (test_angle > angleMax.degrees)
				
					e.soft = false
					e.smooth = false
				end
			end
				
		end
	}
			

end
#reflect points accross z=0 and correct order
def reflect_points(pts_old)
	pts=[]
	
	# reverse order of points
	for i in 0...(pts_old.length)
		pts[i] = pts_old[pts_old.length-1-i]
	end
	
	# reflect z
	for i in 0...(pts.length-1)
		pts[i][2] = -pts[i][2]
	end
	
	return pts
end

model = Sketchup.active_model

model_ents = model.entities

# chamber
p_Rmin = 60.cm
p_Rmax = 160.cm
p_wall_d = 3.cm
p_height = 100.cm

pts_inner = []
pts_inner << [p_Rmin+p_wall_d, 0, p_height-p_wall_d]
pts_inner << [p_Rmax-p_wall_d, 0, p_height-p_wall_d]
pts_inner << [p_Rmax-p_wall_d, 0, -p_height+p_wall_d]
pts_inner << [p_Rmin+p_wall_d, 0, -p_height+p_wall_d]
pts_inner << pts_inner[0]

pts_outer = []

pts_outer << [p_Rmin, 0, p_height]
pts_outer << [p_Rmax, 0, p_height]
pts_outer << [p_Rmax, 0, -p_height]
pts_outer << [p_Rmin, 0, -p_height]
pts_outer << pts_outer[0]


outer_gr = extrude_arc(model_ents, pts_outer, -22.5, 180+22.5, 20)
inner_gr = extrude_arc(model_ents, pts_inner, -22.5, 180+22.5, 20)

chamber_gr=inner_gr.subtract outer_gr

chamber_gr.material='Gray'


# Toroidal Field Coils

p_TF_gap = 5.cm
p_TF_a = 12.cm
p_TF_b = 4.cm
p_TF_bend = 2.cm

h1 = p_height+p_TF_gap
h2 = p_height+p_TF_gap+p_TF_a

r1 = p_Rmax+p_TF_gap
r2 = p_Rmax+p_TF_gap+p_TF_a

pts_outer = []

d = p_TF_gap+p_TF_a

n=20
dT = Math::PI/(2.0*n)
t = 0.0

for i in 0...n
	pts_outer << [p_Rmax + d*Math.cos(t), -p_TF_b, p_height + d*Math.sin(t)]
	t = t + dT
end

for i in 0...n
	pts_outer << [p_Rmin + d*Math.cos(t), -p_TF_b, p_height + d*Math.sin(t)]
	t = t + dT
end

for i in 0...n
	pts_outer << [p_Rmin + d*Math.cos(t), -p_TF_b, -p_height + d*Math.sin(t)]
	t = t + dT
end

for i in 0...n
	pts_outer << [p_Rmax + d*Math.cos(t), -p_TF_b, -p_height + d*Math.sin(t)]
	t = t + dT
end

pts_outer << pts_outer[0]

pts_inner = []

d = p_TF_gap

n=20
dT = Math::PI/(2.0*n)
t = 0.0

for i in 0...n
	pts_inner << [p_Rmax + d*Math.cos(t), -p_TF_b, p_height + d*Math.sin(t)]
	t = t + dT
end

for i in 0...n
	pts_inner << [p_Rmin + d*Math.cos(t), -p_TF_b, p_height + d*Math.sin(t)]
	t = t + dT
end

for i in 0...n
	pts_inner << [p_Rmin + d*Math.cos(t), -p_TF_b, -p_height + d*Math.sin(t)]
	t = t + dT
end

for i in 0...n
	pts_inner << [p_Rmax + d*Math.cos(t), -p_TF_b, -p_height + d*Math.sin(t)]
	t = t + dT
end

pts_inner << pts_inner[0]

outer_gr = model_ents.add_group
outer_ents = outer_gr.entities

outer_face = outer_ents.add_face(outer_ents.add_curve(pts_outer))
outer_face.pushpull -2.0*p_TF_b

inner_gr = model_ents.add_group
inner_ents = inner_gr.entities

inner_face = inner_ents.add_face(inner_ents.add_curve(pts_inner))
inner_face.pushpull -2.0*p_TF_b

tf_gr = inner_gr.subtract outer_gr

soften(tf_gr, 20, 20)

tf_gr.material = 'FireBrick'

n=16
n_max = 9
dT = 2.0*Math::PI/(n)
t = dT

for i in 1...n_max
	next_gr = tf_gr.copy
	next_gr.material = tf_gr.material
	
	trans = Geom::Transformation.rotation [0, 0, 0], [0, 0, 1], t
	
	next_gr.transformation = trans
	
	t = t + dT
end

# Vertical Field Coils
p_VF_gap = 2.cm
p_VF_a = 8.cm
p_VF_b = 4.cm
p_VF_c = 12.cm

r0 = p_Rmax+p_TF_gap+p_TF_a+p_VF_gap

pts = []

pts << [r0, 0, p_VF_b]
pts << [r0+p_VF_a, 0, p_VF_b]
pts << [r0+p_VF_a, 0, -p_VF_b]
pts << [r0, 0, -p_VF_b]
pts << pts[0]

vf_mid_gr = extrude_arc(model_ents, pts, -22.5, 180+22.5, 20)

vf_mid_gr.material = 'CornflowerBlue'

pts = []

pts << [r0, 0, p_height]
pts << [r0+p_VF_a, 0, p_height]
pts << [r0+p_VF_a, 0, p_height-p_VF_c]
pts << [r0, 0, p_height-p_VF_c]
pts << pts[0]

vf_top_gr = extrude_arc(model_ents, pts, -22.5, 180+22.5, 20)

vf_top_gr.material = 'CornflowerBlue'

vf_btm_gr = extrude_arc(model_ents, reflect_points(pts), -22.5, 180+22.5, 20)

vf_btm_gr.material = 'CornflowerBlue'

#viewport

p_vp_a = 40.cm
p_vp_b = 8.cm
p_vp_c = 8.cm
p_vp_d = 0.5.cm

r0 = (p_Rmax + p_Rmin - p_vp_a)/2.0

pts = []
pts << [r0, -p_vp_b/2.0, -p_height]
pts << [r0+p_vp_a, -p_vp_b/2.0, -p_height]
pts << [r0+p_vp_a, p_vp_b/2.0, -p_height]
pts << [r0, p_vp_b/2.0, -p_height]
pts << pts[0]

vp_gr = model_ents.add_group
vp_ents = vp_gr.entities
vp_ents = vp_ents.add_face(vp_ents.add_curve(pts))
vp_ents.pushpull p_wall_d

trans = Geom::Transformation.rotation [0, 0, 0], [0, 0, 1], (180+22.5/2.0).degrees

vp_gr.transformation = trans

chamber_gr = vp_gr.subtract chamber_gr

soften(chamber_gr, 20, 20)
chamber_gr.material='Gray'