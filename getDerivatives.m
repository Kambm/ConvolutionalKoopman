function DB = getDerivatives(B,dt)
    % Gets numerical derivatives of columns of B
    
    DB = zeros(size(B));
    if 1
        % Uses sixth order finite differences
        % for edges derivative
        DB(1:4,:) = (-49/20)*B(1:4,:)+(6)*B(2:5,:)+(-15/2)*B(3:6,:)...
            +(20/3)*B(4:7,:)+(-15/4)*B(5:8,:)+(6/5)*B(6:9,:)+(-1/6)*B(7:10,:);
        DB(end-3:end,:) = (49/20)*B(end-3:end,:)+(-6)*B(end-4:end-1,:)+(15/2)*B(end-5:end-2,:)...
            +(-20/3)*B(end-6:end-3,:)+(15/4)*B(end-7:end-4,:)+(-6/5)*B(end-8:end-5,:)+(1/6)*B(end-9:end-6,:);

        % Uses eighth order finite differences for
        % central derivatives
        DB(5:end-4,:) = (1/280)*B(1:end-8,:)+(-4/105)*B(2:end-7,:)+(1/5)*B(3:end-6,:)...
            +(-4/5)*B(4:end-5,:)+(0)*B(5:end-4,:)+(4/5)*B(6:end-3,:)+(-1/5)*B(7:end-2,:)...
            +(4/105)*B(8:end-1,:)+(-1/280)*B(9:end,:);
    else
        % Uses fifth order finite differences
        % for edges derivative
        DB(1:3,:) = (-137/60)*B(1:3,:)+(5)*B(2:4,:)+(-5)*B(3:5,:)...
            +(10/3)*B(4:6,:)+(-5/4)*B(5:7,:)+(1/5)*B(6:8,:);
        DB(end-2:end,:) = (137/60)*B(end-2:end,:)+(-5)*B(end-3:end-1,:)+(5)*B(end-4:end-2,:)...
            +(-10/3)*B(end-5:end-3,:)+(5/4)*B(end-6:end-4,:)+(-1/5)*B(end-7:end-5,:);
        % Uses sixth order finite differences
        % for central derivatives
        DB(4:end-3,:) = (-1/60)*B(1:end-6,:)+(3/20)*B(2:end-5,:)+(-3/4)*B(3:end-4,:)...
            +(3/4)*B(5:end-2,:)+(-3/20)*B(6:end-1,:)+(1/60)*B(7:end,:);
    end
    DB = DB/dt;
end