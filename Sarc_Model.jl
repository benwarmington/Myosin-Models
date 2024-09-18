
using DifferentialEquations
using GLMakie
using Sundials
using FileIO, JLD2
using ProgressLogging

#GLMakie.activate!();

alpha = 6;
beta = 15;
gamma = 0;
switch = 0;
StS = 0.22;
rot_per = 7.1667*StS;
tooth_per = 6.5*StS;
rot_num = 50;
time = 400;
fil_num = 3;
Alpha = 1.5;
Beta = 0.4*fil_num;
overlap = 1;

# StS = 0.14;
# rot_per = 8.6*StS
# tooth_per = (7.8*StS)/2


tooth_num = Int(ceil((rot_num*rot_per)/tooth_per));
tooth_orig_x=zeros(tooth_num, fil_num);
tooth_loc_x=zeros(tooth_num, fil_num);
tooth_loc_ind=zeros(tooth_num, fil_num);

fil_num*rot_num +1+ fil_num




p = (alpha, beta, switch, StS, rot_per, tooth_per, Int(rot_num), Int(fil_num), Alpha, Beta, time, tooth_num, tooth_orig_x, tooth_loc_x, tooth_loc_ind, overlap, gamma);
u0 = zeros(fil_num*rot_num +1+ fil_num);
du0 = zeros(fil_num*rot_num +1+ fil_num);
N = rot_num;
for k = 1:fil_num
    for i = 1:rot_num

        rot_index = (k-1)*N + k + 1 + i;
        du0[rot_index] = 2*pi;
        u0[rot_index] = 2*pi*rand();
    end
end

display(du0)


function H(theta, StS)

    if mod(theta, 2*pi) > pi - asin(StS) && mod(theta, 2*pi) < pi + asin(StS)

        heavy = 1

    else

        heavy = 0

    end

    return heavy

end


function f2(out, du, u, p, t)

    alpha = p[1]
    beta = p[2]
    StS = p[4]
    rot_per = p[5]
    tooth_per = p[6]
    rot_num = p[7]
    fil_num = p[8]
    Alpha = p[9]
    Beta = p[10]
    Time = p[11]
    tooth_num = p[12]
    tooth_orig_x = p[13]
    tooth_loc_x = p[14]
    tooth_loc_ind = p[15]
    overlap = p[16]
    gamma = p[17]

   
  
    N = rot_num;

   
    
    @inbounds for k = 1:1:fil_num
        @inbounds for j = 1:tooth_num

            tooth_orig_x[j, k] = j*tooth_per - tooth_per;
            tooth_loc_x[j, k] = tooth_orig_x[j, k] + u[(k-1)*N + k + 1];
            tooth_loc_ind[j, k] = mod(tooth_loc_x[j, k], tooth_per*tooth_num);

        end
    end

    rot_forces = zeros(rot_num, fil_num)
    @inbounds for k = 1:1:fil_num
        @inbounds for i = 1:1:rot_num
           
           rot_index = (k-1)*N + k + 1 + i;
        
           atten = (du[rot_index]/(2*pi))*(0.5*tanh(-50*du[rot_index])+0.5);
   
           rot_forces[i, k] = ((1 - du[rot_index]/(2*pi))+atten)*abs(cos(u[rot_index]))*H(u[rot_index], StS);

         end
    end

    springs = zeros(fil_num);

    @inbounds for k = 1:1:fil_num
        
        out[(k-1)*N + k + 1] = sum(@view rot_forces[:, k]) - alpha*(du[(k-1)*N + k + 1]) - beta*(u[(k-1)*N + k + 1] - u[1])

        springs[k] = beta*(u[(k-1)*N + k + 1] - u[1]);
    end

    out[1] = sum(springs)  - Alpha*du[1] - Beta*u[1] - gamma*(0.5*tanh(2*(u[1]-1.3))+0.5)

    @inbounds for k = 1:1:fil_num
        @inbounds for i = 1:1:rot_num
          
            rot_index = (k-1)*N + k + 1 + i;   
            out[rot_index] = 2*pi - du[rot_index]

        end
    end

    @inbounds for k = 1:1:fil_num

        @inbounds for i = 1:rot_num
           
            rot_index = (k-1)*N + k + 1 + i;  

            if H(u[rot_index], StS) == 1 &&  u[(k-1)*N + k + 1] + rot_num*rot_per*overlap > i*rot_per 
                
                @inbounds for j = 1:tooth_num                
                    
                    if abs(tooth_loc_ind[j, k] - ((i-1)*rot_per + sin(-u[rot_index]))) < 0.01 && H(u[rot_index], StS) == 1 

                        out[rot_index] = (du[(k-1)*N + k + 1])/abs(cos(u[rot_index])) - du[rot_index]              

                    end
                
                end

            end
        
        end
    
    end

end


tspan = (0.0, time)
differential_vars = zeros(fil_num*rot_num +1+ fil_num)

@inbounds for i = 1:fil_num*rot_num +1+ fil_num
    differential_vars[i] = true;
end
prob = DAEProblem(f2, du0, u0, tspan, p, differential_vars = differential_vars)

function condition1(out, u, t, integrator) 
 
    rot_per = integrator.p[5]
    tooth_per = integrator.p[6]
    rot_num = integrator.p[7]
    fil_num = integrator.p[8]
    N = rot_num;
   
    tooth_num = integrator.p[12];
    tooth_orig_x = integrator.p[13];
    tooth_loc_x = integrator.p[14];
    tooth_loc_ind = integrator.p[15];
    overlap = integrator.p[16];

    @inbounds for k = 1:1:fil_num
        @inbounds for j = 1:tooth_num

            tooth_orig_x[j, k] = j*tooth_per - tooth_per;
            tooth_loc_x[j, k] = tooth_orig_x[j, k] + u[(k-1)*N + k + 1];
            tooth_loc_ind[j, k] = mod(tooth_loc_x[j, k], tooth_per*tooth_num);

        end
    end
    
    @inbounds for k = 1:fil_num

        @inbounds for i = 1:rot_num
            rot_index = (k-1)*N + k + 1 + i; 
            if H(u[rot_index], integrator.p[4]) == 1 && u[(k-1)*N + k + 1] + rot_num*rot_per*overlap > i*rot_per 

                @inbounds for j = 1:tooth_num                    

                  out[rot_index] = abs(tooth_loc_ind[j, k] - ((i-1)*rot_per + sin(-u[rot_index]))) + 0.01
                
                end
            end
        end
    end

end


function affect!(integrator, idx)
    rot_num = integrator.p[7];
    fil_num = integrator.p[8]
    N = rot_num;

    @inbounds for k = 1:1:fil_num
        @inbounds for i = 1:1:rot_num

            rot_index = (k-1)*N + k + 1 + i;  
            if idx == rot_index
                integrator.du[rot_index] = (integrator.du[(k-1)*N + k + 1])/abs(cos(integrator.u[rot_index]))
            end
        end
    end
end

cb = VectorContinuousCallback(condition1, affect!, fil_num*rot_num +1+ fil_num)

sol = solve(prob, IDA(), maxiters = 10^7, callback=cb, dtmax = 1e-4, reltol = 1e-4, abstol = 1e-6, progress=true, progress_steps=1)


#sol = solve(prob, IDA(linear_solver=:GMRES), maxiters = 20^6, dtmax = 2e-5, reltol = 1e-7, abstol = 1e-8)
#displacement1 = [u[1] for u in sol.u]
#displacement2 = [u[2] for u in sol.u]
#velocity1 = [du[1] for du in sol.u]
#displacement3 = [u[N + 2 + 1] for u in sol.u]
#displacement4 = [u[2*N + 3 + 1] for u in sol.u]
#displacement5 = [u[3*N + 4 + 1] for u in sol.u]
#rot1 = [sin(u[3]) for u in sol.u]
t = sol.t;
# tvals = (time/40000:time/40000:time);
# uvals = sol.(tvals);
# uvals = hcat(uvals...);
displacements = [u[:] for u in sol.u];
displacements = hcat(displacements...);

# vel = zeros(length(tvals), length(uvals[:, 1]));
# bound = zeros(length(tvals), length(uvals[1, :]));

# for j = 1:1:length(uvals[:, 1])
#     for i = 2:1:length(tvals)
#         vel[i, j] = (uvals[j, i] - uvals[j, i-1])/(tvals[i] - tvals[i-1]);
#         if vel[i, j] <  6.2
#             bound[i, j] = 1;
#         end
#     end
# end 
# tot_bound = zeros(length(tvals), 1);
# tot_mot_frac = zeros(length(tvals), 1);
# for i = 1:1:length(tvals)

#     tot_bound[i] = sum(bound[i, :]);
#     tot_mot_frac[i] = tot_bound[i]/(fil_num*rot_num);

# end
# fil_bound = zeros(length(tvals), fil_num);
# fil_mot_frac = zeros(length(tvals), fil_num);
# for j = 1:1:fil_num
#     for i = 1:1:length(tvals)
#         fil_bound[i, j] =  sum(bound[i, (j-1)*rot_num + j + 1:(j-1)*rot_num + j + 1 + rot_num]);
#         fil_mot_frac[i, j] = fil_bound[i, j]/rot_num;
#     end   
# end
# av_vel = zeros(length(tvals))
# window = 60;
# for i = 1:1:length(tvals)-window
#     av_vel[i] = sum(vel[i:i+window, 1])/window;
# end
 
# f = Figure(size = (500, 600))
# ax1 = Axis(f[1, 1], xlabel = "Time", ylabel = "Displacement")
# ax2 = Axis(f[2, 1], xlabel = "Velocity", ylabel = "Force")
# ax3 = Axis(f[3, 1], xlabel = "Time", ylabel = "No. Bound Motors")
# ax4 = Axis(f[4, 1], xlabel = "Time", ylabel = "Sub disp")
# lines!(ax1, t[10:10:end],  displacement1[10:10:end, 1]*(6/StS)) 

# lines!(ax2, vel[10:10:end, 1],  uvals[1, 10:10:end])
# lines!(ax2, av_vel[10:10:end],  uvals[1, 10:10:end])
# lines!(ax3, tvals[10:10:end], tot_mot_frac[10:10:end, 1])
# for j = 1:1:fil_num
#     lines!(ax4, tvals[1:1:end], (uvals[(j-1)*rot_num + j + 1, 1:1:end]-uvals[1, 1:1:end])*(6/StS));
# end
# #lines!(ax4, t[10:10:end], (displacement2[10:10:end, 1] - displacement1[10:10:end, 1])*(6/StS))
# for i = 1:1:fil_num
#     lines!(ax3, tvals[1:1:end], fil_mot_frac[1:1:end, i])
# end
# #lines!(ax3, tvals[2:2:end], fil_mot_frac[2:2:end, 7])

# #lines!(ax1, t[10:10:end], displacement2[10:10:end, 1]*(6/StS))
# #lines!(ax1, t[10:10:end], displacement3[10:10:end, 1]*(6/StS))
# #lines!(ax1, t[10:10:end], displacement4[10:10:end, 1]*(6/StS))
# #lines!(ax1, t[2:2:end], displacement5[2:2:end, 1])
# #lines!(ax1, t[2:2:end],  -rot1[2:2:end, 1])

# display(f)