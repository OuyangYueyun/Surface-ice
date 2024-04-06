function area=caculate_rarea(lat,lon)
    %计算网格占球面面积(r=1)
    n1=size(lat,2);
    n2=size(lon,2);
    area=zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            %纬度范围
            if i==1
                lat1=lat(i);
                lat2=(lat(i)+lat(i+1))/2;
            elseif i==n1
                lat1=(lat(i-1)+lat(i))/2;
                lat2=lat(i);
            else
                lat1=(lat(i)+lat(i-1))/2;
                lat2=(lat(i)+lat(i+1))/2;
            end
            %经度范围
            if j==1
                lon1=lon(j);
                lon2=(lon(j)+lon(j+1))/2;
            elseif j==n2
                lon1=(lon(j-1)+lon(j))/2;
                lon2=lon(j);
            else
                lon1=(lon(j)+lon(j-1))/2;
                lon2=(lon(j)+lon(j+1))/2;
            end
            %面积计算
            area(i,j)=(sin(lat2)-sin(lat1))*(lon2-lon1);
        end
    end
end