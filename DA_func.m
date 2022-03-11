function [gbest,gbestval,fitcount,t]= DA_func(fhd,Dimension,Pop_size,Max_Gen,VRmin,VRmax,X_suru,varargin)
rand('state',sum(100*clock));
Max_iteration=Max_Gen;
SearchAgents_no=Pop_size;
dim=Dimension;
gbest_position=zeros(dim*10000/SearchAgents_no,dim);

ub=VRmax;
lb=VRmin;

if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end

%The initial radius of gragonflies' neighbourhoods
r=(ub-lb)/10;
Delta_max=(ub-lb)/10;

Food_fitness=inf;
Food_pos=zeros(1,dim);

Enemy_fitness=-inf;
Enemy_pos=zeros(1,dim);

X = X_suru;
% X=initialization(SearchAgents_no,dim,ub,lb);
Fitness=zeros(SearchAgents_no,1);
DeltaX=initialization(SearchAgents_no,dim,ub,lb);

% 2. Search History Analysis
if dim==2 % 2 boyut icin yapilacak...
	refresh = 100;
	figure(11)

% initialize
        x1 = linspace(lb(1), ub(1), 101);
        x2 = linspace(lb(2), ub(2), 101);
        x3 = zeros(length(x1), length(x2));
	for i = 1:length(x1)
		for j = 1:length(x2)
			x3(i, j) = feval(fhd,[x1(i);x2(j)],varargin{:});
		end
    end
    
    
% tiles, labels, legend
	str = sprintf('Search History of FN%d',varargin{:});
	xlabel('x_1'); ylabel('x_2');  title(str);
	contour(x1', x2', x3'); hold on;
	plot(X(:,1),X(:,2),'bs','MarkerSize',8);
	drawnow
	% plot(gbest(1),gbest(2),'k*','MarkerSize',8);
end
t=[];iter=1;fitcount=0;
while fitcount<dim*10000 %% && abs(gbestval-varargin{:}*100) > 1e-8
    r=(ub-lb)/4+((ub-lb)*(iter/Max_iteration)*2);
    
    w=0.9-iter*((0.9-0.4)/Max_iteration);
       
    my_c=0.1-iter*((0.1-0)/(Max_iteration/2));
    if my_c<0
        my_c=0;
    end
    
    s=2*rand*my_c; % Seperation weight
    a=2*rand*my_c; % Alignment weight
    c=2*rand*my_c; % Cohesion weight
    f=2*rand;      % Food attraction weight
    e=my_c;        % Enemy distraction weight
    
    for i=1:SearchAgents_no %Calculate all the objective values first
        Fitness(i,1)=feval(fhd,X(i,:)',varargin{:});
		if Fitness(i,1)<Food_fitness
            Food_fitness=Fitness(i,1);
            Food_pos=X(i,:);
        end
        
        if Fitness(i,1)>Enemy_fitness
            if all(X(i,:)<ub) && all( X(i,:)>lb)
                Enemy_fitness=Fitness(i,1);
                Enemy_pos=X(i,:);
            end
        end
    end
    
	fitcount=fitcount+SearchAgents_no;
	
    for i=1:SearchAgents_no
        index=0;
        neighbours_no=0;
        
        clear Neighbours_DeltaX
        clear Neighbours_X
        %find the neighbouring solutions
        for j=1:SearchAgents_no
            Dist2Enemy=distance(X(i,:),X(j,:));
            if (all(Dist2Enemy<=r) && all(Dist2Enemy~=0))
                index=index+1;
                neighbours_no=neighbours_no+1;
                Neighbours_DeltaX(index,:)=DeltaX(j,:);
                Neighbours_X(index,:)=X(j,:);
            end
        end
        
        % Seperation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.1)
        S=zeros(1,dim);
        if neighbours_no>1
            for k=1:neighbours_no
                S=S+(Neighbours_X(k,:)-X(i,:));
            end
            S=-S;
        else
            S=zeros(1,dim);
        end
        
        % Alignment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.2)
        if neighbours_no>1
            A=(sum(Neighbours_DeltaX))/neighbours_no;
        else
            A=DeltaX(i,:);
        end
        
        % Cohesion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.3)
        if neighbours_no>1
            C_temp=(sum(Neighbours_X))/neighbours_no;
        else
            C_temp=X(i,:);
        end
        
        C=C_temp-X(i,:);
        
        % Attraction to food%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.4)
        Dist2Food=distance(X(i,:),Food_pos(1,:));
        if all(Dist2Food<=r)
            F=Food_pos-X(i,:);
        else
            F=zeros(1,dim);
        end
        
        % Distraction from enemy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.5)
        Dist2Enemy=distance(X(i,:),Enemy_pos(1,:));
        if all(Dist2Enemy<=r)
            Enemy=Enemy_pos+X(i,:);
        else
            Enemy=zeros(1,dim);
        end
        
        for tt=1:dim
            if X(i,tt)>ub(tt)
                X(i,tt)=lb(tt);
                DeltaX(i,tt)=rand;
            end
            if X(i,tt)<lb(tt)
                X(i,tt)=ub(tt);
                DeltaX(i,tt)=rand;
            end
        end
        
        if any(Dist2Food>r)
            if neighbours_no>1
                for j=1:dim
                    DeltaX(i,j)=w*DeltaX(i,j)+rand*A(1,j)+rand*C(1,j)+rand*S(1,j);
                    if DeltaX(i,j)>Delta_max(j)
                        DeltaX(i,j)=Delta_max(j);
                    end
                    if DeltaX(i,j)<-Delta_max(j)
                        DeltaX(i,j)=-Delta_max(j);
                    end
                    X(i,j)=X(i,j)+DeltaX(i,j);
                end
            else
                % Eq. (3.8)
                X(i,:)=X(i,:)+levy(dim).*X(i,:);
                DeltaX(i,:)=0;
            end
        else
            for j=1:dim
                % Eq. (3.6)
                DeltaX(i,j)=(a*A(1,j)+c*C(1,j)+s*S(1,j)+f*F(1,j)+e*Enemy(1,j))+w*DeltaX(i,j);
                if DeltaX(i,j)>Delta_max(j)
                    DeltaX(i,j)=Delta_max(j);
                end
                if DeltaX(i,j)<-Delta_max(j)
                    DeltaX(i,j)=-Delta_max(j);
                end
                X(i,j)=X(i,j)+DeltaX(i,j);
            end 
        end
        
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
    end
    gbestval=Food_fitness;
    gbest=Food_pos;   
	
		if fitcount==100*dim||fitcount==200*dim||fitcount==300*dim||fitcount==500*dim...
			||fitcount==1000*dim||fitcount==2000*dim||fitcount==3000*dim||fitcount==4000*dim||fitcount==5000*dim...
			||fitcount==6000*dim||fitcount==7000*dim||fitcount==8000*dim||fitcount==9000*dim||fitcount==10000*dim
				t=[t;abs(gbestval-varargin{:}*100)];              
        end
		gbest_position(iter,:)=gbest;
		Distance = abs(X(1,1)-X(:,1)); % Analiz için eklendi...
		DistanceSum(iter) = sum(Distance)/(SearchAgents_no-1);
		iter=iter+1;
		% 2. Search History Analysis
		if dim==2 && (rem(fitcount/SearchAgents_no,refresh) == 0)% 2 boyut icin yapilacak...
			figure(11)
			plot(X(:,1),X(:,2),'ks','MarkerSize',8,'color',rand(1,3));
		%	plot(gbest(1),gbest(2),'k*','MarkerSize',8);
			xlabel('x_1'); ylabel('x_2'); title(str);
			drawnow
		end
    end
		if dim==2
			% 3. Trajectory Analysis
			figure(12); 
			subplot(2,2,[1 3]);
			contour(x1', x2', x3'); hold on; 
			plot(gbest_position(:,1), gbest_position(:,2), 'kd', 'MarkerSize', 7);
            hold on 
			plot(gbest_position(end,1), gbest_position(end,2), 'rd', 'MarkerSize', 7);
			str = sprintf('Trajectory of Elite FN%d',varargin{:});
			xlabel('x_1'); ylabel('x_2');  title(str);
			subplot(2,2,2);  
			plot(SearchAgents_no:SearchAgents_no:fitcount,gbest_position(:,1),'b');
			xlabel('Function Evaluations'); ylabel('x_1 position');
			subplot(2,2,4);  
			plot(SearchAgents_no:SearchAgents_no:fitcount,gbest_position(:,2),'b');
			xlabel('Function Evaluations'); ylabel('x_2 position');
            drawnow
			% 4. Average Distance Analysis
			figure(13); 
			plot(SearchAgents_no:SearchAgents_no:fitcount,DistanceSum,'b-');
			xlabel('Function Evaluations');
			ylabel('Average distance');
			str = sprintf('Average Distance of FN%d',varargin{:});
            title(str);
            drawnow
		end
end


