using DifferentialEquations
using GLMakie
using Sundials
using ProgressLogging


alpha = 3;
beta = 1;
gamma = 0;
switch = 0;
StS = 0.22;
rot_per = 7.1667*StS;
tooth_per = (6.5*StS)/7;
rot_num = 100;
time = 300;



u0 = zeros(rot_num+1);
du0 = zeros(rot_num+1);

for i = 1:rot_num
    du0[i+1] = 2*pi;
    u0[i+1] = rand()*2*pi;#mod(i*pi, 2*pi);
end





function H(theta, StS)

    if mod(theta, 2*pi) > pi - asin(StS) && mod(theta, 2*pi) < pi + asin(StS)

        heavy = 1

    else

        heavy = 0

    end

    return heavy

end

p = (alpha, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);

function f2(out, du, u, p, t)

    alpha = p[1]
    beta = p[2]
    StS = p[4]
    rot_per = p[5]
    tooth_per = p[6]
    rot_num = Int(p[7])
    gamma = p[8]



    tooth_num = Int(ceil((rot_num*rot_per)/tooth_per));
    tooth_orig_x = zeros(tooth_num);
    tooth_loc_x = zeros(tooth_num);
    tooth_loc_ind = zeros(tooth_num);
    for j = 1:tooth_num
        tooth_orig_x[j] = j*tooth_per - tooth_per;
        tooth_loc_x[j] = tooth_orig_x[j] + u[1];
        tooth_loc_ind[j] = mod(tooth_loc_x[j], tooth_per*tooth_num);
    end


    rot_forces = zeros(rot_num)
  
    for i = 1:1:rot_num
        
        atten =  (du[i+1]/(2*pi))*(0.5*tanh(-80*du[i+1])+0.5);

        rot_forces[i] = ((1 - du[i+1]/(2*pi))+atten)*abs(cos(u[i+1]))*H(u[i+1], StS);

    end


    # if u[1] > 3    
    #     resis = gamma
    #     spring = 0
    # else
    #     resis = 0
    #     spring = beta
    # end


      
        spring = beta


    out[1] = sum(rot_forces)  - alpha*du[1] - spring*u[1] - gamma*(0.5*tanh(2*(u[1]-1.3))+0.5)


    

    for i = 1:rot_num
        out[i+1] = 2*pi - du[i+1]
    end


    for i = 1:rot_num

        if H(u[i+1], StS) == 1


            for j = 7:7:tooth_num


                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.01 && H(u[i+1], 0.1*StS) == 1

                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

                end
            end
            for j = 6:7:tooth_num-1


                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.01 && H(u[i+1], 0.4*StS) == 1

                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

                end
            
            end
            for j = 5:7:tooth_num-2


                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.01 && H(u[i+1], 0.8*StS) == 1

                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

                end
            
            end
            for j = 4:7:tooth_num-3


                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.01 && H(u[i+1], StS) == 1

                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

                end
            
            end
            for j = 3:7:tooth_num-4 

                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.01 && H(u[i+1], 0.8*StS) == 1

                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

                end
            
            end
            for j = 2:7:tooth_num-5

                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.01 && H(u[i+1], 0.4*StS) == 1

                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

                end
            
            end
            for j = 1:7:tooth_num-6

                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.01 && H(u[i+1], 0.1*StS) == 1

                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

                end
            
            end

        
        end
    end



end



function condition1(out, u, t, integrator) 

   
   
    rot_per = integrator.p[5]
    tooth_per = integrator.p[6]
    rot_num = Int(integrator.p[7])

   
    tooth_num = Int(ceil((rot_num*rot_per)/tooth_per));
    tooth_orig_x = zeros(tooth_num);
    tooth_loc_x = zeros(tooth_num);
    tooth_loc_ind = zeros(tooth_num);
    for j = 1:tooth_num
        tooth_orig_x[j] = j*tooth_per - tooth_per;
        tooth_loc_x[j] = tooth_orig_x[j] + u[1];
        tooth_loc_ind[j] = mod(tooth_loc_x[j], tooth_per*tooth_num);
    end

    for i = 1:rot_num

        if H(u[i+1], integrator.p[4]) == 1

            for j = 7:7:tooth_num

                if H(u[i+1], 0.1*integrator.p[4]) == 1        

                    out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.01

                end
            end
            for j = 6:7:tooth_num-1
                
                if H(u[i+1], 0.4*integrator.p[4]) == 1        

                    out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.01

                end

            end
            for j = 5:7:tooth_num-2

                if H(u[i+1], 0.8*integrator.p[4]) == 1        

                    out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.01

                end

            end
            for j = 4:7:tooth_num-3

                if H(u[i+1], integrator.p[4]) == 1        

                    out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.01

                end

            end
            for j = 3:7:tooth_num-4

                if H(u[i+1], 0.8*integrator.p[4]) == 1        

                    out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.01

                end

            end
            for j = 2:7:tooth_num-5

                if H(u[i+1], 0.4*integrator.p[4]) == 1        

                    out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.01

                end

            end
            for j = 1:7:tooth_num-6

                if H(u[i+1], 0.1*integrator.p[4]) == 1        

                    out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.01

                end

            end
        end
    
    end

end


function affect!(integrator, idx)
    rot_num = Int(integrator.p[7]);
    for i = 1:1:rot_num
        if idx == i
            integrator.du[i+1] = integrator.du[1]/abs(cos(integrator.u[i+1]))
        end
    end
end

cb = VectorContinuousCallback(condition1, affect!, rot_num)



tspan = (0.0, time)
differential_vars = zeros(rot_num+1)

for i = 1:rot_num+1
    differential_vars[i] = true;
end

p = (alpha, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);
prob = DAEProblem(f2, du0, u0, tspan, p, differential_vars = differential_vars)

#global p = (alph, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);

sol = solve(prob, IDA(), maxiters = 10^9, callback=cb,  dtmax = 1e-3, reltol = 1e-5, abstol = 1e-6, progress=true, progress_steps=1)
#sol = solve(prob, IDA(linear_solver=:GMRES), maxiters = 20^6, dtmax = 2e-5, reltol = 1e-7, abstol = 1e-8)
displacement = [u[1] for u in sol.u]
rot1 = [sin(u[4]) for u in sol.u]
t = sol.t;
f = Figure(size = (500, 300))
ax1 = Axis(f[1, 1], xlabel = "Time", ylabel = "Displacement")
ax2 = Axis(f[2, 1], xlabel = "Time", ylabel = "No. Bound Motors")
ax3 = Axis(f[3, 1], xlabel = "Time", ylabel = "Velocity")
tvals = (time/20000:time/20000:time);
uvals = sol.(tvals);
uvals = hcat(uvals...);


#displacements = [u[:] for u in sol.u];
#displacements = hcat(displacements...);

vel = zeros(length(tvals), length(uvals[:, 1]));
bound = zeros(length(tvals), length(uvals[1, :]));

for j = 1:1:length(uvals[:, 1])
    for i = 2:1:length(tvals)
        vel[i, j] = (uvals[j, i] - uvals[j, i-1])/(tvals[i] - tvals[i-1]);
        if vel[i, j] <  6.2
            bound[i, j] = 1;
        end
    end
end 
tot_bound = zeros(length(tvals), 1);
tot_mot_frac = zeros(length(tvals), 1);
for i = 1:1:length(tvals)

    tot_bound[i] = sum(bound[i, 2:end]);
    tot_mot_frac[i] = tot_bound[i]/(rot_num);

end
av_vel = zeros(length(tvals))
av_mot_frac = zeros(length(tvals))
window = 1000;
for i = 1:1:length(tvals)-window
    av_vel[i] = sum(vel[i:i+window, 1])/window;
    av_mot_frac[i] = sum(tot_mot_frac[i:i+window])/window;
end
#display(gamma)
display(sum(tot_mot_frac[16001:20000])/4000);
display(sum(vel[16001:20000])/4000); 




lines!(ax1, t[2:2:end],  displacement[2:2:end, 1]*(6/StS)) 
lines!(ax1, t[2:2:end],  -rot1[2:2:end, 1]*(6/StS))
lines!(ax2, tvals[10:10:end], tot_mot_frac[10:10:end, 1])
lines!(ax2, tvals[10:10:end], av_mot_frac[10:10:end, 1])
lines!(ax3, tvals[1:1:end], vel[1:1:end, 1])
lines!(ax3, tvals[1:1:end], av_vel[1:1:end])
#lines!(ax3, uvals[1, 10:10:end], av_vel[10:10:end, 1])
display(f)
