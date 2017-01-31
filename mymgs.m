function V=mymgs(A)
    n=size(A,2);
    V=A;
    
    for i=1:n
        rii=norm(V(:,i));
        qi=V(:,i)/rii;
        for j=i+1:n
            rij=qi'*V(:,j);
            V(:,j)=V(:,j)-rij*qi;
        end
    end
    for j=1:n
        V(:,j)=V(:,j)/norm(V(:,j));
    end
end


    