function value=angleCompute(u,v) 
         value=atan2(norm(cross(u,v)),dot(u,v));
         value=180/pi*value;
end