#import Plots
using GLMakie

function simulate_myosin_system()
    time = 10;
    N =50;
    myo_visc = 0.003;
    fid = 5000*time;
    dev_time = (time/fid):(time/fid):time;
    F_myo = 8;
    fil_visc = 0.00001
    vel_0 = 10000
    F_return = -0.2
    F_spring_const = 0.05;                     
    F_back = 0
    
    h = 0.0000001
 
    F_brownian = 0
    rate_unbind_0 = 70
    rate_bind = 300
    d = 2.75
    p(r, t) = 1 - exp(-r * t)

    act_per = 38.5
    myo_per = 43
    N_act = Int(ceil(((N + 1) * myo_per+2*myo_per) / act_per))
    x_act_orig = zeros(N_act)
    x_act_non_per = similar(x_act_orig)
    

    X0 = zeros(N + 1)
    X0[1] = 0
    for j = 2:N + 1
        X0[j] = (j) * myo_per
    end

    for j = 1:1:N_act
        x_act_orig[j] = j*act_per - act_per;
    end

    
    # for j = 3:3:N_act
    #     x_act_orig[j-2] = j*act_per - 3*act_per-5.5;
    #     x_act_orig[j-1] = j*act_per - 3*act_per;
    #     x_act_orig[j] = j*act_per - 3*act_per + 5.5;

    # end

    t = 0:h:time

    save_step = Int(floor(length(t)/fid));

    x_save = zeros(fid, N + 1)
    vel = zeros(N + 1)
    x_state = ones(N)
    x_state_tracker = ones(Int, fid, N)
    F_elast = zeros(N)
    F_elast_prev = zeros(N) 
    F_elast_save = zeros(fid, N)     
    x_head = zeros(N+1);
    bind_site = zeros(N);
    F_brown = zeros(fid, N)
    F_const = 0
    x_act = zeros(N_act)
    x_previous = zeros(N+1)
    x_current = zeros(N+1)
    x_next = zeros(N+1)
    x_previous .= X0[:]
    x_current .= X0[:]
    x_head[1:N] .= X0[2:N+1]
    #x_head[2, 1:N] .= X0[2:N+1]

    for c = 2:length(t) - 1
        F_BROWNIAN = zeros(N)
        F_MYO = zeros(N)
        F_RETURN = zeros(N)
        F_HOLD = zeros(N)
        L_diff = zeros(N)
        MYO_VISC = zeros(N)

        if isapprox(mod(t[c], 1), 0.0)
            println(t[c])
        end

        length(t)/fid

        if mod(c, save_step) == 0.0

           i = Int(c/save_step);
           x_save[i , :] = x_current;
           x_state_tracker[i, :] = x_state;
           F_elast_save[i, :] = F_elast;

        end

        if c == 20000
            F_const = F_back
        end

        for k = 1:N_act
            x_act_non_per[k] = x_act_orig[k] + x_current[1]
            x_act[k] = mod(x_act_non_per[k], (act_per * N_act)-39)
        end

        for j = 1:N
            if x_state[j] == 1
               x_head[j] = X0[j+1]
                for k = 1:N_act
                    if abs(x_current[j + 1] - x_act[k]) < d
                        P_binding = p(rate_bind, h)*exp(-(abs(x_current[j+1] - x_act[k]))^2/(2*(d/3)^2))
                        if P_binding > rand()
                            #x[c, j + 1] = x_act[k]
                            bind_site[j] = k
                            x_state[j] = 2
                        end
                    end
                end
            elseif x_state[j] == 2
                rate_unbind = 0
                x_head[j] = x_act[Int(bind_site[j])]
                rate_unbind = rate_unbind_0 * (0.907 * exp(-abs(F_elast_prev[j]) * 0.362) + 0.093 * exp(0.121 * abs(F_elast_prev[j])))
                if rand() < p(rate_unbind, h) || abs(X0[j + 1] - x_head[j]) > 2 * myo_per
                    x_state[j] = 4
                elseif x_current[j + 1] >= X0[j + 1] + 8

        
                    x_state[j] = 3
                end
            elseif x_state[j] == 3
                rate_unbind = 0
                x_head[j] = x_act[Int(bind_site[j])]
                rate_unbind = rate_unbind_0 * (1.257 * exp(-abs(F_elast_prev[j]) * 0.6036) + 0.107 * exp(0.0966 * abs(F_elast_prev[j])))
                if rand() < p(rate_unbind, h) || abs(x_current[j + 1] - x_head[j]) > 2 * myo_per
                    x_state[j] = 4
                end
            elseif x_state[j] == 4
                x_head[j] = x_current[j + 1]
                if x_current[j + 1] <= X0[j + 1]
                    x_state[j] = 1
                end
            end
        end

        #vel[N + 1] .= (x_next[N + 1] - x_current[N + 1]) / h

        for j = 1:N
            L_diff[j] = (x_current[j + 1] - x_head[j])
            F_elast[j] = (0.03579 * L_diff[j] + 1.732 * exp(1.3765 * L_diff[j]) - 1.535)
            if x_state[j] == 1
                F_MYO[j] = 0
                F_RETURN[j] = 0
                F_BROWNIAN[j] = F_brownian * randn()
                F_HOLD[j] = 0
                MYO_VISC[j] = 0
                #x_state_tracker[c, j] = 1
            elseif x_state[j] == 2
                F_MYO[j] = F_myo
                F_RETURN[j] = 0
                F_BROWNIAN[j] = F_brownian * randn()
                F_HOLD[j] = 0
                MYO_VISC[j] = myo_visc
               # x_state_tracker[c, j] = 2
            elseif x_state[j] == 3
                F_MYO[j] = 0
                F_RETURN[j] = 0
                F_BROWNIAN[j] = F_brownian * randn()
                F_HOLD[j] = F_elast[j]
                MYO_VISC[j] = myo_visc
              #  x_state_tracker[c, j] = 3
            elseif x_state[j] == 4
                F_MYO[j] = 0
                F_RETURN[j] = F_return
                F_BROWNIAN[j] = F_brownian * randn()
                F_HOLD[j] = 0
                MYO_VISC[j] = 0
              #  x_state_tracker[c, j] = 4
            end
        end

        motor_force = zeros(N)
        for j = 1:N
            if x_state[j] == 2 || x_state[j] == 3
                motor_force[j] = copy(F_elast[j])
            end
        end

        x_next[1] = (sum(motor_force) + F_spring_const * X0[1] + ((4 * fil_visc * x_current[1]) / (2 * h)) - ((fil_visc * x_previous[1]) / (2 * h)) - F_const) / ((3 * fil_visc) / (2 * h) + F_spring_const)
        vel[1] = (x_next[1] - x_current[1]) / h

        F_all = zeros(N)
        for j = 1:N
            F_all[j] = - F_elast[j] + F_BROWNIAN[j] + F_RETURN[j] + F_MYO[j] + F_HOLD[j]
            if vel[1] > 0
                 x_next[j + 1] = (1 / 3) * (4 * x_current[j + 1] - x_previous[j + 1] + 2 * h * (F_all[j] / ((F_MYO[j] / vel_0) + fil_visc + MYO_VISC[j])))
            else
                x_next[j + 1] = (1 / 3) * (4 * x_current[j + 1] - x_previous[j + 1] + 2 * h * (F_all[j] / (fil_visc + MYO_VISC[j]))) 
            end
        end
        x_previous, x_current, x_next = x_current, x_next, x_previous;
        F_elast_prev = F_elast;
    end
    return x_save, x_state_tracker, F_elast_save, dev_time, N, X0, F_spring_const

end


#  for N = 10:5:100
# # # Call the function
#      track = zeros(5)
#      for i = 1:1:5
#         x, x_state_tracker, vel, F_elast, t, x_act = simulate_system(N)

#         println("N = $(N) test $(i)   -   $(x[end, 1]/1.5)")
#         track[i] = x[end, 1]/1.5
#      end
    
#      println("N= $(N) Average        $(sum(track)/5)")
    
#  end
# Plotting

x, x_state_tracker, F_elast, t,  N, X0, F_spring_const = simulate_myosin_system()
#p1 = CairoMakie.lines(t, x[:, 5:9], lw = 2)
#p2 = CairoMakie.lines(t, F_elast[:, 9:9], lw = 2)
GLMakie.activate!()
p3 = lines(t, x[:, 1], lw = 3)

bound = zeros(length(t), length(x[1, :])-1)
for i = 1:1:length(t)
    for j = 1:1:length(x[1, :])-1
         if (x_state_tracker[i, j] == 2 || x_state_tracker[i, j] == 3) && F_elast[i, j] > 0
            bound[i, j] = 1
         end
    end
end


f2 = Figure(size = (300, 350))
ax7 = Axis(f2[1, 1], xlabel = "Time (s)", ylabel = "Myosin Number")
p4 = heatmap!(ax7, bound[3005:5:4000, :])
save("Late_Wave.png", f2, px_per_unit = 2) 
# tot_bound = zeros(length(t));

# N = length(x[1, :])-1

# for i = 1:1:length(t)
#     bound = 0;
#     for j = 1:1:N
#          if x_state_tracker[i, j] == 2 || x_state_tracker[i, j] == 3
#             bound = bound+1;
#          end
#     end
#     tot_bound[i] = bound;
# end

# p4 = Plots.plot(t[200:200:end], (tot_bound[200:200:end]/(N))*100, lw = 2)

# vel_smooth = zeros(length(vel[:, 1]),1);
# w = Int(8000);
# for i = 1:length(vel[:, 1])-(w+10)
#     s = sum(vel[i:i+w, 1]);
#     vel_smooth[i] = (s/(w+1));
#     # if vel_smooth[i] < 0
#     #     vel_smooth[i] = 0;
#     # end

# end
# p5 = Plots.plot(t[1:end-w], vel_smooth[1:end-w, 1]);


# Display the plots
 #display(p1)
# display(p2)
 display(p3)
#display(f2)

#display(p4)
#display(p5)