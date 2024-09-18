# These packages allow the use of the required solvers - the SUNDIALS suite of solver packages is particularly useful for solving non-linear and discontinuous equations.
# Makie is a better visualisation package than Plots, the base julia data visualisation package.

using DifferentialEquations
using GLMakie
using Sundials
using ProgressLogging
using FFTW
GLMakie.activate!();

alpha = 0.2;
beta = 0;
gamma = 0.2;
switch = 0;
# StS = 0.8;
#StS = 0.22;
#StS = 0.18
StS = 0.9;
# rot_per = 7.1667*StS;
# tooth_per = (6.5*StS);

#Bound for L_t checking
#rot_per = 1.7;#8.6*StS
#tooth_per = 0.1;#(7.8*StS)

rot_per = 8.6*StS
tooth_per = (7.8*StS)/3

rot_per = 3
tooth_per = 2.2


#StS_Checking param!
#rot_per = 1.85
#tooth_per = (1.75)


rot_num = 11;
time = 100;

# if rot_per > tooth_per
#     phase_diff = 2*pi - (mod(rot_per, tooth_per)/tooth_per)*2*pi;
# elseif tooth_per > rot_per
#     phase_diff = (mod(tooth_per, rot_per)/tooth_per)*2*pi;
# else
#     phase_diff = 0;
# end


u0 = zeros(rot_num+1);
du0 = zeros(rot_num+1);

for i = 1:rot_num
    du0[i+1] = 2*pi;
    u0[i+1] = 0;#rand()*2*pi#mod(phase_diff*i, 2*pi);
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
    @inbounds for j = 1:tooth_num
        tooth_orig_x[j] = j*tooth_per - tooth_per;
        tooth_loc_x[j] = tooth_orig_x[j] + u[1];
        tooth_loc_ind[j] = mod(tooth_loc_x[j], tooth_per*tooth_num);
    end


    rot_forces = zeros(rot_num)
  
    # @inbounds for i = 1:1:rot_num
        
    #     atten =  (du[i+1]/(2*pi))*(0.5*tanh(-80*du[i+1])+0.5);

    #     rot_forces[i] = ((1 - du[i+1]/(2*pi))+atten)*abs(cos(u[i+1]))*H(u[i+1], StS);

    # end

    
    for i = 1:1:rot_num

        velad = du[i+1]-1.5*pi;
        
        atten =  (velad/(0.5*pi))*(0.5*tanh(-80*velad)+0.5);

        rot_forces[i] = ((1 - velad/(0.5*pi))+atten)*abs(cos(u[i+1]))*H(u[i+1], StS);

    end


    # if u[1] > 300    
    #      spring = beta;
        
    # # elseif t > 50 && t < 100
    # #     resis = 0

    # else 

    #     spring = 0;

    # end


      
       spring = beta


    out[1] = sum(rot_forces)  - alpha*du[1] - spring*(u[1])- gamma*(0.5*tanh(2*(u[1]-1.3))+0.5)


    

    @inbounds for i = 1:rot_num
        out[i+1] = 2*pi - du[i+1]
    end


    @inbounds for i = 1:rot_num

        @inbounds for j = 1:tooth_num


            if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.02 && H(u[i+1], StS) == 1

                out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              

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
    @inbounds for j = 1:tooth_num
        tooth_orig_x[j] = j*tooth_per - tooth_per;
        tooth_loc_x[j] = tooth_orig_x[j] + u[1];
        tooth_loc_ind[j] = mod(tooth_loc_x[j], tooth_per*tooth_num);
    end

    @inbounds for i = 1:rot_num

        if H(u[i+1], integrator.p[4]) == 1        

            @inbounds for j = 1:tooth_num

                out[i] = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) + 0.02

            end
        end
    end

end


function affect!(integrator, idx)
    rot_num = Int(integrator.p[7]);
    @inbounds for i = 1:1:rot_num
        if idx == i
            integrator.du[i+1] = integrator.du[1]/abs(cos(integrator.u[i+1]))
        end
    end
end

cb = VectorContinuousCallback(condition1, affect!, rot_num)


Force_Vel_Data = zeros(20, 3)

#for c = 4:1:20
c = 1;
    #StS = c*0.02;
    # alpha = 2;

    # rot_per = 1
   # tooth_per = 0.11*c

   # alpha = c;
    
    if rot_per > tooth_per
        phase_diff =  mod(((tooth_per - rot_per)/tooth_per)*2*pi, 2*pi);
    elseif tooth_per > rot_per
        phase_diff =   mod(((tooth_per -  rot_per)/tooth_per)*2*pi, 2*pi);
    else
        phase_diff = 0;
    end
    
    #phase_diff = ((2*pi)/40)*c;

    u0 = zeros(rot_num+1);
    du0 = zeros(rot_num+1);
    
    for i = 1:rot_num
        du0[i+1] = 2*pi;
        u0[i+1] = 0;#mod(i*phase_diff, 2*pi);
    end
    






   # alpha = gamma; 
    tspan = (0.0, time)
    differential_vars = zeros(rot_num+1)

    @inbounds for i = 1:rot_num+1
        differential_vars[i] = true;
    end

    p = (alpha, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);
    prob = DAEProblem(f2, du0, u0, tspan, p, differential_vars = differential_vars)

    #global p = (alph, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);

    sol = solve(prob, IDA(), maxiters = 10^9, callback=cb,  dtmax = 1e-4, reltol = 1e-4, abstol = 1e-6, progress=true, progress_steps=1)
    #sol = solve(prob, IDA(linear_solver=:GMRES), maxiters = 20^6, dtmax = 2e-5, reltol = 1e-7, abstol = 1e-8)
  
    #rot1 = [cos(u[4]) for u in sol.u]
    #t = sol.t;
    f = Figure(size = (500, 300))
    ax1 = Axis(f[1, 1], xlabel = "Time", ylabel = "Displacement")
    ax2 = Axis(f[2, 1], xlabel = "Time", ylabel = "No. Bound Motors")
    ax3 = Axis(f[3, 1], xlabel = "Time", ylabel = "Velocity")
    ax4 = Axis(f[4, 1], xlabel = "Time", ylabel = "x pos")
    tvals = (time/20000:time/20000:time);
    uvals = sol.(tvals);
    uvals = hcat(uvals...);


    #displacements = [u[:] for u in sol.u];
    #displacements = hcat(displacements...);

    vel = zeros(length(tvals), length(uvals[:, 1]));
    bound = zeros(length(tvals), length(uvals[1, :]));

    @inbounds for j = 1:1:length(uvals[:, 1])
        @inbounds for i = 2:1:length(tvals)
            vel[i, j] = (uvals[j, i] - uvals[j, i-1])/(tvals[i] - tvals[i-1]);
            if vel[i, j] <  6.2
                bound[i, j] = 1;
            end
        end
    end 
    tot_bound = zeros(length(tvals), 1);
    tot_mot_frac = zeros(length(tvals), 1);
    @inbounds for i = 1:1:length(tvals)

        tot_bound[i] = sum(bound[i, 2:end]);
        tot_mot_frac[i] = tot_bound[i]/(rot_num);

    end
    av_vel = zeros(length(tvals))
    av_mot_frac = zeros(length(tvals))
    window = 100;
    for i = 1:1:length(tvals)-window
        av_vel[i] = sum(vel[i:i+window, 1])/window;
        av_mot_frac[i] = sum(tot_mot_frac[i:i+window])/window;
    end
    display(gamma)
    display(sum(tot_mot_frac[16001:20000])/4000);
    display(sum(vel[16001:20000])/4000); 
   phase_diff = zeros(rot_num-1, 1)
    for i = 1:1:rot_num-1
        phase_diff[i] = mod(uvals[i+2, 19900] - uvals[i+1, 19900], 2*pi)
    end
    av_phase_diff = sum(phase_diff[1:end])/length(phase_diff[1:end])
    display(av_phase_diff)


    
     index = c;
    
    
   Force_Vel_Data[index, 3] = sum(vel[16001:20000])/4000;
   Force_Vel_Data[index, 2] = sum(tot_mot_frac[16001:20000])/4000;
   Force_Vel_Data[index, 1] = alpha;



    n = length(uvals[1, :])
    fs = (20000/time)
    frequencies = (0:n-1)*(fs/n)



    #sig_mean = sum(Disps[1, :])/length(Disps[1, :])

    fft_result = fft(cos.(uvals[30, :] ))
    magnitude = abs.(fft_result)


    normalised_magnitude = magnitude/n

    lines!(ax2, frequencies[1:800], normalised_magnitude[1:800], color = :blue, linewidth = 2)

    (ind, loc) = findmax(normalised_magnitude[1:800])
    p_freq = frequencies[loc]
    phase_diff =  ((tooth_per - rot_per)/tooth_per)*2*pi;
    wave_per = ((2*pi)/abs(phase_diff))*rot_per

    wave_prop_vel = p_freq*wave_per
    display(wave_prop_vel)
#end



#FileIO.save("C:/Users/yp20984/OneDrive - University of Bristol/Julia/TestJulia/Vel/LTsweep_LR_2_75.jld2", "Data", Force_Vel_Data)

# for gamma = 5.25:0.25:25

#     time = 100;

#     alpha = gamma; 
#     tspan = (0.0, time)
#     differential_vars = zeros(rot_num+1)

#     for i = 1:rot_num+1
#         differential_vars[i] = true;
#     end

#     p = (alpha, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);
#     prob = DAEProblem(f2, du0, u0, tspan, p, differential_vars = differential_vars)

#     #global p = (alph, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);

#     sol = solve(prob, IDA(), maxiters = 10^9, callback=cb,  dtmax = 5e-4, reltol = 1e-6, abstol = 1e-7, progress=true, progress_steps=1)
#     #sol = solve(prob, IDA(linear_solver=:GMRES), maxiters = 20^6, dtmax = 2e-5, reltol = 1e-7, abstol = 1e-8)
#     displacement = [u[1] for u in sol.u]
#     rot1 = [sin(u[4]) for u in sol.u]
#     t = sol.t;
#     f = Figure(size = (500, 300))
#     ax1 = Axis(f[1, 1], xlabel = "Time", ylabel = "Displacement")
#     ax2 = Axis(f[2, 1], xlabel = "Time", ylabel = "No. Bound Motors")
#     ax3 = Axis(f[3, 1], xlabel = "Time", ylabel = "Velocity")
#     tvals = (time/20000:time/20000:time);
#     uvals = sol.(tvals);
#     uvals = hcat(uvals...);


#     #displacements = [u[:] for u in sol.u];
#     #displacements = hcat(displacements...);

#     vel = zeros(length(tvals), length(uvals[:, 1]));
#     bound = zeros(length(tvals), length(uvals[1, :]));

#     for j = 1:1:length(uvals[:, 1])
#         for i = 2:1:length(tvals)
#             vel[i, j] = (uvals[j, i] - uvals[j, i-1])/(tvals[i] - tvals[i-1]);
#             if vel[i, j] <  6.2
#                 bound[i, j] = 1;
#             end
#         end
#     end 
#     tot_bound = zeros(length(tvals), 1);
#     tot_mot_frac = zeros(length(tvals), 1);
#     for i = 1:1:length(tvals)

#         tot_bound[i] = sum(bound[i, 2:end]);
#         tot_mot_frac[i] = tot_bound[i]/(rot_num);

#     end
#     av_vel = zeros(length(tvals))
#     av_mot_frac = zeros(length(tvals))
#     window = 1000;
#     for i = 1:1:length(tvals)-window
#         av_vel[i] = sum(vel[i:i+window, 1])/window;
#         av_mot_frac[i] = sum(tot_mot_frac[i:i+window])/window;
#     end
#     display(gamma)
#     display(sum(tot_mot_frac[16001:20000])/4000);
#     display(sum(vel[16001:20000])/4000); 


#     index = Int(gamma*4)
    
#    Force_Vel_Data[index, 3] = sum(vel[16001:20000])/4000;
#    Force_Vel_Data[index, 2] = sum(tot_mot_frac[16001:20000])/4000;
#    Force_Vel_Data[index, 1] = gamma;
# end



lines!(ax1, tvals[2:2:end], uvals[1, 2:2:end]*beta) 
#lines!(ax1, t[2:2:end],  displacement[2:2:end, 1]*(0.015) .+ gamma) 
lines!(ax1, tvals[2:2:end],  cos.(uvals[2, 2:2:end]))
lines!(ax2, tvals[10:10:end], tot_mot_frac[10:10:end, 1])
lines!(ax2, tvals[10:10:end], av_mot_frac[10:10:end, 1])
lines!(ax3, tvals[1:1:end], vel[1:1:end, 1])
lines!(ax3, tvals[1:1:end], av_vel[1:1:end])
#lines!(ax4, 1:50, sin.(uvals[2:51, 16000]))
#heatmap!(tvals[2000:2:20000], rot_per.*(1:50), cos.(transpose(uvals[2:51, 2000:2:20000])))
#lines!(ax3, uvals[1, 10:10:end], av_vel[10:10:end, 1])

display(f)

#FileIO.save("/user/work/yp20984/Vel/Vel_for_RL_2_25.jld2", "Data", Force_Vel_Data)

#display(av_mot_frac[end-window-5]);
#display(av_vel[end-window-5]);
# fig2 = Figure(size = (180, 130))
# ax2 = Axis(fig2[1, 1], xlabel = "L_t", ylabel = "No. Bound Motors")
# lines!(ax2, 0.11*5:0.11:20*0.11, Force_Vel_Data[5:20, 2], linewidth = 2, color = :black)
# display(fig2)
# save("Cont_L_t.png", fig2, px_per_unit = 8)