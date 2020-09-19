function [S] = orthogonal_reflection(V, P, n)
%% orthogonal_reflection : function to compute the -orthogonal- orthogonal_reflectionion of a vector
% or an array of vectors toward an hyperplane of the 3D or the 2D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2020.
%
%
% Syntax
%
% S = orthogonal_reflection(V, P, n);
%
%
% Description
%
% S = orthogonal_reflection(V, P, n) computes the vector S, which results from
% the orthogonal orthogonal_reflectionion of vector V toward the hyperplane P.
%
%
% Input arguments
%
%       [ -Vx- ]
% - V = [ -Vy- ], real (array of) column vector(s) double, the vector(s) to orthogonal_reflection. Size(V) = [3,vector_nb]. 
%       [ -Vz- ] 
%
%         [Px]
% - P : = [Py], real column vector double, one point belonging to the orthogonal_reflectionion hyperplane. Size(P) = [3,1].
%         [Pz]
%
%         [nx]
% - n : = [ny], real column vector double, one vector normal to the orthogonal_reflectionion hyperplane. Size(n) = [3,1].
%         [nz]
%
%
% Output argument
%
%       [ -Sx- ]
% - S = [ -Sy- ], real (array of) vector(s) double, the resulting orthogonal_reflectioned vector(s). Size(S) = [3,vector_nb].
%       [ -Sz- ]
%
%
% Example #1 : 2D point
% n = [1 1]';
% P = n;
% V = [2 2];
% S = orthogonal_reflection(V, P, n) % expected : S = [0; 0]
%
%
% Example #2 : array of 2D points
% n = [2 -1]';
% P = [1 0]';
% nb_pts = 16;
% V = rand(2,nb_pts) - cat(1,zeros(1,nb_pts),ones(1,nb_pts));
% S = orthogonal_reflection(V,P,n);
%
% figure;
% plot(V(1,:),V(2,:),'b+','MarkerSize',8,'Linewidth',2), hold on;
% plot(S(1,:),S(2,:),'g+','MarkerSize',8,'Linewidth',2), hold on;
% line([0 1],[-2 0],'Color',[1 0 0],'Linewidth',2), hold on;
% axis equal, axis tight;
% set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);
%
%
% Example #3 : array of 3D points 1
% n = [1 1 1]';
% P = [0 0 0]';
% nb_pts = 32;
% V = rand(3,nb_pts);
% S = orthogonal_reflection(V,P,n);
%
% figure;
% [X,Y] = meshgrid(-1:0.1:1,1:-0.1:-1);
% Z = -X-Y;
% surf(X,Y,Z), shading interp, hold on;
%
% plot3(V(1,:),V(2,:),V(3,:),'b+','MarkerSize',8,'Linewidth',2), hold on;
% plot3(S(1,:),S(2,:),S(3,:),'g+','MarkerSize',8,'Linewidth',2), hold on;
% cellfun(@(c1,c2) line([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'LineStyle','--','Color',[1 1 1]),num2cell(V,1),num2cell(S,1),'un',0);
%
% colormap([1 0 0]);
% alpha(0.5);
% axis equal, axis tight;
% view(45,26);
% set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);
%
%
% Example #4 : array of 3D points 2
% n = [1 1 1]';
% P = [0 0 0]';
% V = cat(1,[2 1.25 1 1 1.25 1.75 2 2 1.75 1],zeros(1,10),[2 2 1.75 1.25 1 1 0.75 0.25 0 0]);
% S = orthogonal_reflection(V,P,n);
%
% figure;
% [X,Y] = meshgrid(-1:0.1:1,1:-0.1:-1);
% Z = -X-Y;
% surf(X,Y,Z), shading interp, hold on;
%
% plot3(V(1,:),V(2,:),V(3,:),'b+','MarkerSize',8,'Linewidth',2), hold on;
% line(V(1,:),V(2,:),V(3,:),'Color',[0 0 1],'Linewidth',2), hold on;
% plot3(S(1,:),S(2,:),S(3,:),'g+','MarkerSize',8,'Linewidth',2), hold on;
% line(S(1,:),S(2,:),S(3,:),'Color',[0 1 0],'Linewidth',2), hold on;
%
% colormap([1 0 0]);
% alpha(0.5);
% axis square, axis tight;
% view(42,35);
% set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);


%% Input parsing
assert(nargin > 2, 'Not enough input arguments.');
assert(nargin < 4, 'Too many input arguments.');

% Format check
assert(isnumeric(V) && isreal(V), 'Inputs V must be real (array of) column vector(s).');
assert(isnumeric(n) && isreal(n) && size(n,2) == 1 && norm(n) > 0, 'Inputs n must be real column non null vector.');

Ndim = size(V,1);
nb_pts = size(V,2);
assert(Ndim == 2 || Ndim == 3, 'Dimensions accepted must be an integer in the range |[2, 3]|.');


%% Body
switch Ndim 
    
    case 2
        
        u = cat(2,-n(2,1),n(1,1)); % one line director vector
        [~,H] = point_to_line_distance([0 0],u,P');
        
    case 3
       
        [~,H] = point_to_plane_distance([0 0 0],n',P');       
        
end

n = n / sqrt(sum(n.^2,1));
Sm = eye(Ndim) - 2*(n*n');
S = Sm * V + 2*repmat(H',[1,nb_pts]);
S(abs(S) < 1e4*eps) = 0;


end % orthogonal_reflection


%% Point to line distance subfunction
function [d2H, H] = point_to_line_distance(P, u, I0)


nb_pts = size(P,1);
dimension = size(P,2);    

% Zeros padding in 2D case
if dimension == 2
        
    t_H = (u(1)*(P(:,1)-repmat(I0(1),[nb_pts,1])) + ...
           u(2)*(P(:,2)-repmat(I0(2),[nb_pts,1]))) / ...
           sum(u.^2);
    
    x_H = I0(1) + t_H*u(1);
    y_H = I0(2) + t_H*u(2);    
    
    % Orthogonal projected point
    H = zeros(size(P,1),dimension);
    H(:,1) = x_H;
    H(:,2) = y_H;    
    
    % Distance
    d2H = sqrt(sum((P-H).^2,2));
    H = H(:,1:dimension);
    
elseif dimension == 3
    
    t_H = (u(1)*(P(:,1)-repmat(I0(1),[nb_pts,1])) + ...
           u(2)*(P(:,2)-repmat(I0(2),[nb_pts,1])) + ...
           u(3)*(P(:,3)-repmat(I0(3),[nb_pts,1])) ) / ...
           sum(u.^2);
    
    x_H = I0(1) + t_H*u(1);
    y_H = I0(2) + t_H*u(2);
    z_H = I0(3) + t_H*u(3);
    
    % Orthogonal projected point
    H = zeros(size(P,1),dimension);
    H(:,1) = x_H;
    H(:,2) = y_H;
    H(:,3) = z_H;
    
    % Distance
    d2H = sqrt(sum((P-H).^2,2));
    H = H(:,1:dimension);
    
end


end % point_to_line_distance


%% Point to plane distance subfunction
function [d2H, H] = point_to_plane_distance(M, n, I)


nb_pts = size(M,1);

d_I = -(n(1)*I(1)+n(2)*I(2)+n(3)*I(3));
t_H = -(repmat(d_I,[nb_pts,1])+n(1)*M(:,1)+n(2)*M(:,2)+n(3)*M(:,3)) / sum(n.^2);

x_H = M(:,1) + t_H*n(1);
y_H = M(:,2) + t_H*n(2);
z_H = M(:,3) + t_H*n(3);

H = cat(2,x_H,y_H,z_H);
d2H = sqrt(sum((M-H).^2,2));


end % point_to_plane_distance