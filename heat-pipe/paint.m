function paint(A,dx,dr)
contourf(A)
grid on
grid minor
xticks([1 1/dx])
xStart = string(dx/2);
xEnd = string(1-dx/2);
xticklabels({xStart,xEnd})
yticks([1 (0.8-0.2)/dr])
rStart = string(0.2 + dr/2);
rEnd = string(0.8 - dr/2);
yticklabels({rStart,rEnd})
end