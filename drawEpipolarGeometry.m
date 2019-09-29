function drawEpipolarGeometry(p1,p2,E,K,imageSize)
    run('/home/hongyi/Documents/Package/vlfeat/vlfeat/toolbox/vl_setup.m');
    F = transpose(inv(K)) * E * inv(K);
    
    x = 1:imageSize(2);
    for i = 1:size(p1,2)
        l2 = F * p1(:,i);
        y2(i,:) = (-l2(1) .* x - l2(3)) / l2(2);
    end
    
    for i = 1:size(p2,2)
        l1 = F' * p2(:,i);
        y1(i,:) = (-l1(1) .* x - l1(3)) / l1(2);
    end
    vl_tightsubplot(1,2,1,'spacing',0.05);
    colormap = rand(5,3);
    for i = 1:size(p1,2)
        plot(p1(1,:),p1(2,:),'k+','MarkerSize',15); hold on
        plot(x, y1(i,:), '-','color',colormap(i,:)); hold on
    end
    hold off
    xlim([1 imageSize(2)]);
    ylim([1 imageSize(1)]);
    %axis equal
    hold off
    vl_tightsubplot(1,2,2,'spacing',0.05);
    for i = 1:size(p2,2)
        plot(p2(1,:),p2(2,:),'k+','MarkerSize',15); hold on
        plot(x, y2(i,:), '-','color',colormap(i,:)); hold on
    end
    xlim([1 imageSize(2)]);
    ylim([1 imageSize(1)]);
    %axis equal
    hold off
end