clear; clc;
load('SensorData10.mat');  data = SensorData10;
dt = mean(diff(data(:,1)));

time = data(:,1);
ref_att = data(:,2:4);      % 参考姿态 [pitch, roll, yaw]
ref_vel = data(:,5:7);      % 参考速度 [VE, VN, VU]
ref_pos = data(:,8:10);     % 参考位置 [Lat, Lon, H]

omega_ie = 7.2921151467e-5;  R_e = 6378137.0;
f = 1/298.257223563;          e2 = 2*f - f^2;
g0 = 9.7803267714;

gyro = data(:,11:13);  acc = data(:,14:16);
N = size(data,1);
n_sub = 2;  % 双子样

% ===== 初始姿态四元数: q = q_Z(yaw) * q_X(pitch) * q_Y(roll) =====
att0 = deg2rad(data(1,2:4));
p = att0(1); r = att0(2); y = att0(3);

cr = cos(r/2); sr = sin(r/2);
cp = cos(p/2); sp = sin(p/2);
cy = cos(y/2); sy = sin(y/2);

a = [cp; sp; 0; 0];  b = [cr; 0; sr; 0];
q_xy = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
        a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
        a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
        a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)];

a = [cy; 0; 0; sy];  b = q_xy;
q = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
     a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
     a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
     a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)];
q = q / norm(q);

vel = data(1,5:7)';
lat = deg2rad(data(1,8));  lon = deg2rad(data(1,9));  h = data(1,10);

att = zeros(N,3);  vel_out = zeros(N,3);  pos_out = zeros(N,3);
att(1,:) = rad2deg([p, r, y]);
vel_out(1,:) = vel';
pos_out(1,:) = [rad2deg(lat), rad2deg(lon), h];

% ===== 捷联惯导解算主循环 =====
for k = 1:N-1
    % 导航参数
    sinL = sin(lat);  cosL = cos(lat);
    R_M = R_e*(1-e2)/(1-e2*sinL^2)^1.5 + h;
    R_N = R_e/sqrt(1-e2*sinL^2) + h;
    w_ie_n = [0; omega_ie*cosL; omega_ie*sinL];
    w_en_n = [-vel(2)/R_M; vel(1)/R_N; vel(1)*tan(lat)/R_N];
    w_in_n = w_ie_n + w_en_n;

    % ===== 步长起始时刻四元数和姿态矩阵（保存备用） =====
    q_old = q;
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    Cbn_old = [q0^2+q1^2-q2^2-q3^2,  2*(q1*q2-q0*q3),      2*(q0*q2+q1*q3);
               2*(q1*q2+q0*q3),       q0^2-q1^2+q2^2-q3^2,  2*(q2*q3-q0*q1);
               2*(q1*q3-q0*q2),       2*(q0*q1+q2*q3),      q0^2-q1^2-q2^2+q3^2];

    w_in_b = Cbn_old' * w_in_n;

    % ===== 双子样旋转矢量 + 圆锥误差补偿 =====
    dtheta_sub = zeros(3, 2);
    for i = 1:2
        alpha = (i - 0.5) / 2;
        w_ib_b_interp = (1 - alpha) * gyro(k,:)' + alpha * gyro(k+1,:)';
        dtheta_sub(:,i) = (w_ib_b_interp - w_in_b) * (dt / 2);
    end
    dtheta = dtheta_sub(:,1) + dtheta_sub(:,2) ...
           + (2/3) * cross(dtheta_sub(:,1), dtheta_sub(:,2));

    % 四元数更新
    phi = norm(dtheta);
    if phi > 1e-12
        dq = [cos(phi/2); dtheta/phi*sin(phi/2)];
    else
        dq = [1; dtheta/2];  dq = dq/norm(dq);
    end

    a = q;  b = dq;
    q = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
         a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
         a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
         a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)];
    q = q / norm(q);

    % ===== 姿态提取 =====
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    Cbn = [q0^2+q1^2-q2^2-q3^2,  2*(q1*q2-q0*q3),      2*(q0*q2+q1*q3);
           2*(q1*q2+q0*q3),       q0^2-q1^2+q2^2-q3^2,  2*(q2*q3-q0*q1);
           2*(q1*q3-q0*q2),       2*(q0*q1+q2*q3),      q0^2-q1^2-q2^2+q3^2];

    C32 = Cbn(3,2);  C31 = Cbn(3,1);  C33 = Cbn(3,3);
    C12 = Cbn(1,2);  C22 = Cbn(2,2);
    pitch = asin(C32);
    roll  = atan2(-C31, C33);
    yaw   = -atan2(-C12, C22);
    if roll < 0, roll = roll + 2*pi; end
    att(k+1,:) = rad2deg([pitch, roll, yaw]);

    % ===== 双子样速度增量 =====
    dtheta_v1 = gyro(k,:)'   * (dt/2);
    dtheta_v2 = gyro(k+1,:)' * (dt/2);

    dv_sub = zeros(3, 2);
    for i = 1:2
        alpha = (i - 0.5) / 2;
        acc_interp = (1 - alpha) * acc(k,:)' + alpha * acc(k+1,:)';
        dv_sub(:,i) = acc_interp * (dt / 2);
    end

    % 标准双子样划桨补偿
    dv_rot  = 0.5 * cross(dtheta_v1 + dtheta_v2, dv_sub(:,1) + dv_sub(:,2));
    dv_scul = 0.5 * (cross(dtheta_v1, dv_sub(:,2)) + cross(dv_sub(:,1), dtheta_v2));
    dv_b = dv_sub(:,1) + dv_sub(:,2) + dv_rot + dv_scul;

    % ===== 中间时刻姿态矩阵：从q_old旋转半步 =====
    dtheta_half = (dtheta_v1 + dtheta_v2) / 2;
    phi_h = norm(dtheta_half);
    if phi_h > 1e-12
        dq_h = [cos(phi_h/2); dtheta_half/phi_h*sin(phi_h/2)];
    else
        dq_h = [1; dtheta_half/2];  dq_h = dq_h/norm(dq_h);
    end
    a = q_old;  b = dq_h;
    q_mid = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
             a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
             a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
             a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)];
    q_mid = q_mid / norm(q_mid);

    q0=q_mid(1); q1=q_mid(2); q2=q_mid(3); q3=q_mid(4);
    Cbn_mid = [q0^2+q1^2-q2^2-q3^2,  2*(q1*q2-q0*q3),      2*(q0*q2+q1*q3);
               2*(q1*q2+q0*q3),       q0^2-q1^2+q2^2-q3^2,  2*(q2*q3-q0*q1);
               2*(q1*q3-q0*q2),       2*(q0*q1+q2*q3),      q0^2-q1^2-q2^2+q3^2];

    % 比力用中间时刻姿态矩阵投影
    dvel_n = Cbn_mid * dv_b;

    g_h = g0*(1+0.00527094*sinL^2+2.32718e-5*sinL^4) - 3.086e-6*h;
    g_n = [0; 0; -g_h];
    coriolis = cross(2*w_ie_n + w_en_n, vel);

    vel_new = vel + dvel_n + (g_n - coriolis) * dt;

    % ===== 位置更新（梯形积分，经度分母各用各时刻值） =====
    lat_new = lat + 0.5 * (vel(2) + vel_new(2)) / R_M * dt;
    sinL_new = sin(lat_new);  cosL_new = cos(lat_new);
    R_N_new  = R_e/sqrt(1-e2*sinL_new^2) + h;
    lon = lon + 0.5 * (vel(1)/(R_N*cosL) + vel_new(1)/(R_N_new*cosL_new)) * dt;
    lat = lat_new;
    h   = h + 0.5 * (vel(3) + vel_new(3)) * dt;

    vel = vel_new;
    vel_out(k+1,:) = vel';
    pos_out(k+1,:) = [rad2deg(lat), rad2deg(lon), h];
end

% ===== 结果输出 =====
fprintf('最终姿态: Pitch=%.2f°, Roll=%.1f°, Yaw=%.1f°\n', att(end,:));
fprintf('最终速度: VE=%.2f, VN=%.2f, VU=%.2f m/s\n', vel_out(end,:));
fprintf('最终位置: Lat=%.6f°, Lon=%.6f°, H=%.2f m\n', pos_out(end,:));

fprintf('\n===== RMS 误差 =====\n');
att_err = att - ref_att;
att_err(:,2) = mod(att_err(:,2)+180, 360) - 180;
att_err(:,3) = mod(att_err(:,3)+180, 360) - 180;
fprintf('姿态: Pitch=%.4f°, Roll=%.4f°, Yaw=%.4f°\n', rms(att_err));
vel_err = vel_out - ref_vel;
fprintf('速度: VE=%.4f, VN=%.4f, VU=%.4f m/s\n', rms(vel_err));
pos_err = pos_out - ref_pos;
fprintf('位置: Lat=%.6f°, Lon=%.6f°, H=%.4f m\n', rms(pos_err));

% ===== 绘图 =====
plotData = {att(:,1), ref_att(:,1); att(:,2), ref_att(:,2); att(:,3), ref_att(:,3);
            vel_out(:,1), ref_vel(:,1); vel_out(:,2), ref_vel(:,2); vel_out(:,3), ref_vel(:,3);
            pos_out(:,1), ref_pos(:,1); pos_out(:,2), ref_pos(:,2); pos_out(:,3), ref_pos(:,3)};
ylabels  = {'Pitch (°)','Roll (°)','Yaw (°)', ...
            'VE (m/s)','VN (m/s)','VU (m/s)', ...
            '纬度 (°)','经度 (°)','高度 (m)'};
titles   = {'俯仰角','滚转角','偏航角', ...
            '东向速度','北向速度','天向速度', ...
            '纬度','经度','高度'};

for i = 1:9
    figure('Name', titles{i});
    plot(time, plotData{i,1}, 'b', time, plotData{i,2}, 'r--', 'LineWidth',1.2);
    ylabel(ylabels{i}); xlabel('时间 (s)');
    legend('惯导','参考'); grid on; title(titles{i});
end