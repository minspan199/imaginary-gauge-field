%plot ring-like intensity distribbution figure
function [F]=edgeplot1(x0,y0,z,f)
t=50;    % 50*50 pixels per ring
x=(x0-1)*t+(1:t); y=(y0-1)*t+(1:t);


for i=(x0-1)*t+(1:t)  
     for j=(y0-1)*t+(1:t)
          if   ((i-t/2-(x0-1)*t)^2+(j-t/2-(y0-1)*t)^2)>10^2      % ring inner radius 10 pixels
              if ((i-t/2-(x0-1)*t)^2+(j-t/2-(y0-1)*t)^2)<15^2     % ring outter radius 15 pixels
                  f(1+abs(round(i)),1+abs(round(j)))=z;
              end
          end
     end
end
  F=f;        

