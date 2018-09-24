function color_line_plot(x, colorscale)
    if size(x,2) < size(x,1)
        x = x';
    end
    if size(colorscale,2) < size(colorscale,1)
        colorscale = colorscale';
    end
    colorscale = colorscale - min(colorscale);
    max(colorscale)
    colorscale = colorscale/max(colorscale);
    negcolor = [0,0,1];
    poscolor = [0,1,0];
    
    figure;
    hold on;
    if size(x,1) == 1
        t = 1:size(x,2);
        for j = 2:size(x,2)
            j
            plot(t(j-1:j),x(j-1:j), 'color',colorscale(j)*poscolor+(1-colorscale(j))*negcolor, 'linewidth', 2);
        end
    elseif size(x,1) == 2
        for j = 2:size(x,2)
            j
            plot(x(1,j-1:j),x(2,j-1:j),'color',colorscale(j)*poscolor+(1-colorscale(j))*negcolor, 'linewidth', 2);
        end
    elseif size(x,1) == 3
        for j = 2:size(x,2)
            j
            plot3(x(1,j-1:j),x(2,j-1:j),x(3,j-1:j),'color',colorscale(j)*poscolor+(1-colorscale(j))*negcolor, 'linewidth', 2);
        end
    end
end