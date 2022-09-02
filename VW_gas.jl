#=This software will run a molecular dynamics simulation for atoms with a LJ potential. 
The variables are defined in the following way:
l is the length of the box
w is the width of the box
h is the height of the box
sigma is defined for LJ potential
epsilon is defined for LJ potential
time is the total simulation steps
dt is the time increment=#

using LinearAlgebra
boxsides = 16
l = 10
w = 10
h = 10
num_particles = l*w*h
dim = [l,w,h]
coord = zeros(Float64,l*w*h ,3)
prev_coord = zeros(Float64,l*w*h ,3)
forces = zeros(Float64,l*w*h,3)
dt = 0.01
sigma = 1.0;
epsilon = 1.0;
temperature = 1.0
#println(coord)
#initialCoords(coord,1,dim)
#print(coord)
# here we start the output file. The information included will be the number of molecules and the name of the molecules. We will then use
# append to add on the trajectory frames. 
outfile = "trajectory.xyz"
f = open(outfile,"w")
println(f,string(num_particles))
println(f,"LJ_gas")
close(f)



#Initialize the starting coordinates for the gas on an equal lattice spacing.
#currently it only allows a grid spacing of 1.0 but will be updated to allow for any spacing. 
function initialCoords!(coordinates, spacing, dim)
    length = dim[1]
    width = dim[2]
    height  = dim[3]
    
    for i in 1:length
        for j in 1:width
            for k in 1:height
                coordinates[length*height*(i-1) + height*(j-1) + k,:] = [i, j, k] 
            end
        end
    end
end


function initVelocities!(temperature,coordinates,prev_coordinates,time_step)
        num_elements = div(length(coordinates),3)
        #I need to check this degree of freedome calc.
        nFreedom = 3*num_elements-3
        velocities = rand(-0.5:0.0002:0.5,1000,3)
        #prev_coordinates = velocities
        #remove the center of mass motion of the system
        sumV = [0.0,0.0,0.0]
        for i in 1:3
            sumV[i] = sum(velocities[:,i])/num_elements
        end
        for i in 1:num_elements
            velocities[i,:] -= sumV
        end
        #calculate the KE and rescale the velocities to match the temperature based on an gas
        v2 = 0.0
        for i in 1:num_elements
            v2 += norm(velocities[i,:])^2
        end
        scale = sqrt(nFreedom*temperature/v2)
        velocities*=scale
        for i in 1:num_elements
            prev_coordinates[i,:] = coordinates[i,:] - velocities[i,:]*time_step
        end
end
#Define the lennard jones energy
function energyLJ(r,sigma,epsilon)
    4.0*epsilon*((sigma/r)^6-(sigma/r)^12)
end
#Define the lennard jones energy shifted to allow for a cut-off distance
function energyLJ_shift(r, rc,sigma,epsilon)
        if r<rc
            energyLJ(r,sigma,epsilon)-energyLJ(rc,sigma,epsilon)
        else
            0.0
	end
end

#Define the function to compute the forces between the particles
function computeForces!(forces,coordinates, boxlength, rc)
    fill!(forces,0.0)  
    num_elements = div(length(coordinates),3)
    #num_elements = 1000
    for i in 1:num_elements - 1
        for j in i+1:num_elements
            distVec = coordinates[i,:] - coordinates[j,:]
            #apply minimum image convention
            distVec[1] -=boxlength*round(distVec[1]/boxlength)
            distVec[2] -=boxlength*round(distVec[2]/boxlength)
            distVec[3] -=boxlength*round(distVec[3]/boxlength)
            dist = norm(distVec)
            #only need to calculate a force if it is within the cut-off distance
            if dist <= rc
                pairwiseforce = (48.0*dist^(-14))-(24.0*dist^(-8))
                forces[i,:] += pairwiseforce*distVec
                forces[j,:] -= pairwiseforce*distVec
            end
        end
    end
end

#using the verlet algorithm to integrate the EOM
function integrateEOM!(coordinates,prev_coordinates,forces,time_step,boxLength)
    temp_coordinates = [0.0, 0.0, 0.0] 
    num_elements = div(length(coordinates),3)
    for i in 1:num_elements
        temp_coordinates = (2.0*coordinates[i,:] - prev_coordinates[i,:] + (time_step^2)*forces[i,:]).%boxLength
        prev_coordinates[i,:] = coordinates[i,:]
        coordinates[i,:]= temp_coordinates
    end
end
#initialCoords(coord,1,dim)
#print(coord)
#println(energyLJ(0.5,sigma,epsilon))
#println(energyLJ_shift(5.1,3,sigma,epsilon))
#println(forces)
#println(coord[1,:]-coord[2,:])
#println(forces)
#velocities = rand(-.5:0.0001:0.5,1000,3)
#prev_coord = velocities
#println(prev_coord)
#println(coord)


#outfile = "trajectory.xyz"
#initialize simulation
initialCoords!(coord,1,dim)
initVelocities!(temperature,coord,prev_coord,dt)
for m in 1:100
    computeForces!(forces,coord,l,2.5)
    integrateEOM!(coord,prev_coord,forces,dt,boxsides)
    f = open(outfile,"a")
    for i in 1:num_particles
        println(f,"Atom"*string(i)*" "*string(coord[i,1])*" "*string(coord[i,2])*" "*string(coord[i,3]))
    end
    close(f) 
end
#close(f)
