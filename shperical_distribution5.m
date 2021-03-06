clear
%% initialization
K = 1;
fi = @(theta1,phi1,theta2,phi2) K / ((cos(phi1)*sin(theta1)-cos(phi2)*sin(theta2))^2+(sin(phi1)*sin(theta1)-sin(phi2)*sin(theta2))^2+(cos(theta1)-cos(theta2))^2)^0.5;
fitheta = @(theta1,phi1,theta2,phi2) K / ((cos(phi1)*sin(theta1)-cos(phi2)*sin(theta2))^2+(sin(phi1)*sin(theta1)-sin(phi2)*sin(theta2))^2+(cos(theta1)-cos(theta2))^2)^1.5 * ((cos(phi1)*sin(theta1)-cos(phi2)*sin(theta2))*cos(phi1)*cos(theta1)+(sin(phi1)*sin(theta1)-sin(phi2)*sin(theta2))*sin(phi1)*cos(theta1)-(cos(theta1)-cos(theta2))*sin(theta1));
fiphi = @(theta1,phi1,theta2,phi2) K / ((cos(phi1)*sin(theta1)-cos(phi2)*sin(theta2))^2+(sin(phi1)*sin(theta1)-sin(phi2)*sin(theta2))^2+(cos(theta1)-cos(theta2))^2)^1.5 * (-(cos(phi1)*sin(theta1)-cos(phi2)*sin(theta2))*sin(phi1)*sin(theta1)+(sin(phi1)*sin(theta1)-sin(phi2)*sin(theta2))*cos(phi1)*sin(theta1));

N_MAX = 80;
mean_area = zeros(N_MAX,1);
max_area = zeros(N_MAX,1);
min_area = zeros(N_MAX,1);
for N=4:N_MAX
    N
    rng(23 ,'twister');
    theta = unifrnd(0,pi,[N,1]);%rand(N,1)*pi;
    phi = unifrnd(0,2*pi,[N,1]);%rand(N,1)*2*pi;
    MAXITERS = 3000;
    dt = 0.01;
    
    Px0=sin(theta).*cos(phi);
    Py0=sin(theta).*sin(phi);
    Pz0=cos(theta);
    
    
    %% main loop
    for iter = 1:MAXITERS
        xstar = zeros(2,N);
        for i=1:N
            x0 = [theta(i);phi(i)];
            dx_theta = 0;
            dx_phi = 0;
            for j=1:N
               if (i==j)
                   continue;
               end
               dx_theta = fitheta(x0(1),x0(2),theta(j),phi(j)) + dx_theta;
               dx_phi = fiphi(x0(1),x0(2),theta(j),phi(j)) + dx_phi;
               dx = [dx_theta;dx_phi];
            end
            xstar(:,i) = x0 + dx*dt;
        end
        x = xstar;
        for i=1:N
            if (x(1,i)<0)
                x(1,i) = -x(1,i);
                x(2,i) = x(2,i) + pi;
            elseif (x(1,i)>pi)
                x(1,i) = 2*pi - x(1,i);
                x(2,i) = x(2,i) + pi;
            end
            if (x(2,i)>2*pi)
                x(2,i) = x(2,i) - 2*pi;
            elseif (x(2,i)<0)
                x(2,i) = x(2,i) + 2*pi;
            end
        end
        theta = x(1,:)';
        phi = x(2,:)';
    end
    
    %% discussion
    Px=sin(theta).*cos(phi);
    Py=sin(theta).*sin(phi);
    Pz=cos(theta);
    [K,Vol]=convhull(Px,Py,Pz);
    %scatter3(Px,Py,Pz,'filled')
    area = zeros(size(K,1),1);
    for i = 1:size(K,1)
        v1 = [Px(K(i,1))-Px(K(i,2));Py(K(i,1))-Py(K(i,2));Pz(K(i,1))-Pz(K(i,2))];
        v2 = [Px(K(i,1))-Px(K(i,3));Py(K(i,1))-Py(K(i,3));Pz(K(i,1))-Pz(K(i,3))];
        area(i) = 0.5*sqrt(norm(v1)^2*norm(v2)^2-sum(v1.*v2)^2);
    end
    
    mean_area(N) = mean(area);
    %std(area)
    max_area(N) = max(area);
    min_area(N) = min(area);
    
end

max_area = max_area(4:N_MAX,1);
min_area = min_area(4:N_MAX,1);
mean_area = mean_area(4:N_MAX,1);
%mean_area_function = @(n) sqrt(3)/2*(1-2*cos(n*pi/(3*n-6)))/(1-cos(n*pi/(3*n-6)));

%% drawing result
%figure
%scatter3(Px0,Py0,Pz0,'filled')
%hold on
%trisurf(K,Px,Py,Pz);
%hold on
%sphere
