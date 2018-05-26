function [errors] = project(width,height,v_top,v_bottom,v_left,v_right,x_0,y_0,vx_0,vy_0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

flag=0;
h=20*height;
w=20*width;
b=height;
a=width;
q_m=-1.76*10^11;

%Sets up original matrix, with wall potentials set
bcs=[v_bottom,v_top,v_right,v_left];
potential=zeros(h,w);
potential(1,:)=v_top;
potential(h,:)=v_bottom;
potential(:,1)=v_left;
potential(:,w)=v_right;

%creates original residual matrix      
residual=zeros(h,w);

%loop runs until tolerance has been met in all element of residual matrix
while(flag == 0)

    %updates residual matrix based on current potential matrix
    for i=2:(h-1)
        for j=2:(w-1)
            residual(i,j)=potential(i-1,j)+potential(i+1,j)+potential(i,j-1)+potential(i,j+1)-4*potential(i,j);
        end
    end

    %updates potential matrix based on current residual matrix
    for i=2:(h-1)
        for j=2:(w-1)
            potential(i,j)=potential(i,j)+0.25*residual(i,j);
            %original(i,j)=round(original(i,j));
        end
    end
    flag=test(residual);
end

%function tests to see if tolerance has been met in residual matrix
%returns 1 if has been satisfied, 0 if not
function [result] = test(A)
    result=1;
    for i=2:(h-1)
        for j=2:(w-1)
            if(abs(A(i,j))>0.0002)
                result=0;
            end
        end
    end
end

%Makes heatmap of potential
subplot(3,1,1), imagesc(potential), colorbar;


%Makes 25 term approximation for potential 
analytic_potential_step=zeros(h,w);
A=zeros(1,50);
B=zeros(1,50);
C=zeros(1,50);
D=zeros(1,50);
for n=1:2:50
    A(n)=(4/(n*pi))*((bcs(3)-bcs(4)*exp(-n*pi*a/b))/(1-exp(-2*n*pi*a/b)));
    B(n)=(4*exp(-n*pi*a/b)/(n*pi))*((bcs(4)-bcs(3)*exp(-n*pi*a/b))/(1-exp(-2*n*pi*a/b)));
    C(n)=(4/(n*pi))*((bcs(2)-bcs(1)*exp(-n*pi*b/a))/(1-exp(-2*n*pi*b/a)));
    D(n)=(4*exp(-n*pi*b/a)/(n*pi))*((bcs(1)-bcs(2)*exp(-n*pi*b/a))/(1-exp(-2*n*pi*b/a)));
    for i=2:(h-1)
        for j=2:(w-1)
            x=width-(j/20);
            y=i/20;
            horizontal=((A(n)*exp(-n*pi*x/b)+B(n)*exp(n*pi*x/b))*sin(n*pi*y/b));
            vertical=((C(n)*exp(-n*pi*y/a)+D(n)*exp(n*pi*y/a))*sin(n*pi*x/a));
            analytic_potential_step(i,j)=(vertical)+(horizontal);
        end
    end
    
    if n==1
        analytic_potential=analytic_potential_step;
    else
        analytic_potential=analytic_potential+analytic_potential_step;
    end
    
end

%Uses coefficients from above to form analytic E-field for point (x,y)
function [E] = Efield(x,y)
    Ef_x=0;
    Ef_y=0;
    for n=1:2:5
        Ef_x=Ef_x+((A(n)*exp(-n*pi*x/b)-B(n)*exp(n*pi*x/b))*(n*pi/b)*sin(n*pi*y/b))-(n*pi/a)*((C(n)*exp(-n*pi*y/a)+D(n)*exp(n*pi*y/a))*cos(n*pi*x/a));
        Ef_y=Ef_y-((A(n)*exp(-n*pi*x/b)+B(n)*exp(n*pi*x/b))*cos(n*pi*y/b))*(n*pi/b)+(n*pi/a)*((C(n)*exp(-n*pi*y/a)-D(n)*exp(n*pi*y/a))*sin(n*pi*x/a));
        E=[Ef_x, Ef_y];
    end
end

%calculate error between relaxed and anayltic matrices
errors = error_calc(analytic_potential, potential);

%graphs on heatmap the analytic potential
analytic_potential(1,:)=v_top;
analytic_potential(h,:)=v_bottom;
analytic_potential(:,1)=v_left;
analytic_potential(:,w)=v_right;
subplot(3,1,2), imagesc(analytic_potential), colorbar;


%tracking of particle in field (pos=[x,y] vel=[vx,vy] acc=[ax,ay])
t_f=10^-4;
t_step=10^-8;
i=1;
t(1)=0;
time=t(1);
pos(1,:)=[width-x_0,height-y_0];
vel(1,:)=[-vx_0,-vy_0];

while ( (bound(x,y) == 0) && (time <= t_f))
   
    i = i + 1;
    time = time + t_step;
    t(i)=time;
 
    kpos1=pos(i-1,:);
    kvel1=vel(i-1,:);
    kacc1=(q_m)*Efield(pos(i-1,1),pos(i-1,2));
    
    kpos2=pos(i-1,:)+0.5*kvel1*t_step;
    kvel2=vel(i-1,:)+0.5*kacc1*t_step;
    kacc2=(q_m)*Efield(kpos2(1),kpos2(2));
    
    kpos3=pos(i-1,:)+0.5*kvel2*t_step;
    kvel3=vel(i-1,:)+0.5*kacc2*t_step;
    kacc3=(q_m)*Efield(kpos3(1),kpos3(2));
    
    kpos4=pos(i-1,:)+kvel3*t_step;
    kvel4=vel(i-1,:)+kacc3*t_step;
    kacc4=(q_m)*Efield(kpos4(1),kpos4(2));
    
    pos(i,:)=pos(i-1,:)+(t_step/6)*(kvel1+2*kvel2+2*kvel3+kvel4);
    vel(i,:)=vel(i-1,:)+(t_step/6)*(kacc1+2*kacc2+2*kacc3+kacc4);
    
end
subplot(3,1,3), plot((width-pos(:,1)),(height-pos(:,2)));
subplot(3,1,3), axis([0 width 0 height]);


%function for calculating error (calculates summed total and rms error)
function [error] = error_calc(A,B)
    err=0;
    for i=2:(h-1)
        for j=2:(w-1)
           err=err+(A(i,j)-B(i,j))^2;
        end
    end
    rms=sqrt(err/(i*j));
    err_tot=err;
    error=[rms,err_tot];
end

%function determining if the particle has hit a wall or not
function [hit] = bound(x,y)
    hit=0;
    if((x<0)||(x>width))
        hit=1;
    elseif((y<0)||(y>height))
        hit=1;
    end
end

end

