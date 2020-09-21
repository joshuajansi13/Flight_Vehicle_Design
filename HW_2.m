clc
clear all

rho = 1.225;
v_infinity = 3:0.5:20;
wing_span = 0.9144;
tail_span = 0.1524;
wing_chord = 0.20352;
tail_chord = 0.1016;
wing_thickness = 0.0508;
tail_thickness = 0.0381;

reference_area = wing_span*wing_chord;
K_constant = 0.38;
aspect_ratio = wing_span^2 / reference_area;

fuselage_diameter = 0.1;
fuselage_length = 0.6096;
wing_thickness_to_chord_ratio = wing_thickness/wing_chord;
tail_thickness_to_chord_ratio = tail_thickness/tail_chord;

sweep_of_wing = 0;
sweep_of_tail = 0;

e_inv = 0.98*(1-2*(fuselage_diameter/wing_span)^2);

mu = 1.81*10^-5;

exposed_area_of_wing = wing_span*wing_chord;
exposed_area_of_tail = tail_span*tail_chord;

wetted_areas_of_wing_tail_fuselage = [
    2*(1+0.2*wing_thickness_to_chord_ratio)*exposed_area_of_wing;
    2*(1+0.2*tail_thickness_to_chord_ratio)*exposed_area_of_tail;
    0.75*pi*fuselage_diameter*fuselage_length
    ];

k_values_of_wing_tail_fuselage = [
    1+2*cos(sweep_of_wing)*wing_thickness_to_chord_ratio+100*wing_thickness_to_chord_ratio^4;
    1+2*cos(sweep_of_tail)*tail_thickness_to_chord_ratio+100*tail_thickness_to_chord_ratio^4;
    calculate_k_of_fuselage(fuselage_length, fuselage_diameter)
    ];

lift_of_plane_in_newtons = 9.32;

for index = 1:length(v_infinity)
    dynamic_pressure = 0.5*rho*v_infinity(index)^2;
    
    reynolds_numbers_wing_tail_fuselage = [
        rho*v_infinity(index)*wing_chord/mu;
        rho*v_infinity(index)*tail_chord/mu;
        rho*v_infinity(index)*fuselage_length/mu
        ];
    
    C_f_turbulent_wing_tail_fuselage = [
        0.074/reynolds_numbers_wing_tail_fuselage(1)^0.2;
        0.074/reynolds_numbers_wing_tail_fuselage(2)^0.2;
        0.074/reynolds_numbers_wing_tail_fuselage(3)^0.2
        ];
    
    C_d_p_of_wing_tail_fuselage = [
        k_values_of_wing_tail_fuselage(1)*C_f_turbulent_wing_tail_fuselage(1)*wetted_areas_of_wing_tail_fuselage(1)/reference_area;
        k_values_of_wing_tail_fuselage(2)*C_f_turbulent_wing_tail_fuselage(2)*wetted_areas_of_wing_tail_fuselage(2)/reference_area;
        k_values_of_wing_tail_fuselage(3)*C_f_turbulent_wing_tail_fuselage(3)*wetted_areas_of_wing_tail_fuselage(3)/reference_area
        ];
    
    total_parasitic_drag_coefficient = C_d_p_of_wing_tail_fuselage(1)+C_d_p_of_wing_tail_fuselage(2)+C_d_p_of_wing_tail_fuselage(3);
    
    e = 1 / ((1/e_inv)+K_constant*total_parasitic_drag_coefficient*pi*aspect_ratio);
   
    total_drag(index) = total_parasitic_drag_coefficient*dynamic_pressure*reference_area+lift_of_plane_in_newtons^2 / (dynamic_pressure*pi*wing_span^2*e);
    
    lift_to_drag_ratio(index) = lift_of_plane_in_newtons/total_drag(index);
end

[lift_to_drag_ratio_max, index] = max(lift_to_drag_ratio);
corresponding_v_infinity = v_infinity(index);

figure(1)
plot(v_infinity, lift_to_drag_ratio, 'b-')
hold on
plot(corresponding_v_infinity, lift_to_drag_ratio_max, 'kx', 'MarkerSize', 10)
xlabel('Velocity (m/s)')
ylabel('Lift-to-Drag Ratio')

% TODO: find design flight speed (speed that maximizes L/D)- v_infinity_max
%       corresponding Re
%       corresponding lift coefficient - corresponding_lift_coefficient
%       max L/D - lift_to_drag_ratio_max

corresponding_dynamic_pressure = 0.5*rho*v_infinity(index)^2;
corresponding_lift_coefficient = lift_of_plane_in_newtons/(corresponding_dynamic_pressure*reference_area);

function [k_fuse] = calculate_k_of_fuselage(fuselage_length, fuselage_diameter)
    fineness_ratio = fuselage_length/fuselage_diameter;
    if fineness_ratio < 15 && fineness_ratio > 5
        k_fuse = 1.675-0.09*fineness_ratio+0.003*fineness_ratio^2;
    elseif fineness_ratio >= 15
        k_fuse = 1;
    else
        fprintf("Error");
    end
end

