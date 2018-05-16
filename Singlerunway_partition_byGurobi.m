clear all
clc

format

total_part=0; % 最後總共拆成幾組
gurobi_time=0; % gurobi的計算時間
res_type=[]; % 飛機經過Algo之後的類別序列
res_num_prev=[]; % 飛機經過Algo之後各自最後的進場順序
res_num=[]; % 最後經過調換後，進場的編號順序
num_tmp = 0; % 進場順序的tmp紀錄器

%% 初始值
n_f=8;    % Number of aircraft
%  X= randint(1,n_f,[1 3]) ;
X=[1 1 3 1 3 3 2 2];
 E_T =sort(ceil((n_f*10)*10*rand(1,n_f)));
% E_T=[248,264,265,299,366,771,823,845,856,1065,1128,1288,1365,1443,1475];
EST=[E_T-E_T(1)]';
 v=ceil(50*rand(1,n_f))+480;  %巡航速度 %節去計算
% v=[505,507,490,512,482,496,507,497,511,499,487,526,513,513,514];
 D=ceil(28*rand(1,n_f))+150  ; %降落速度
% D=[171,174,165,155,161,173,155,160,178,173,157,178,152,162,162];

Res=zeros(n_f); %拆成n_f組
partition= 1; %紀錄總共有幾組
tmp_head = 1; %紀錄每一組的第一台飛機是幾號飛機
ptr = 2; %暫時的指標

t0 = clock; % 開始計時
%% 目的在將 1122232123 可以拆解成-> 1 12 2 232123 這些類別 ， 每組數量越少algorithm總計算時間就會下降

for i=2:length(X)
    if ptr == length(X)
        if  X(ptr) ~= X(ptr-1)
            Xlen = ptr-tmp_head+1;
            Res(partition , 1:Xlen) = X(tmp_head:ptr);
        elseif X(ptr) == X(ptr-1)
            if tmp_head == ptr-1
                Res(partition, 1:2) = X(ptr-1:ptr);
                break;
            end
            Xlen = ptr-tmp_head;
            Res(partition , 1:Xlen) = X(tmp_head:(ptr-1));
            partition = partition + 1;
            Res(partition , 1) = X(ptr);
        end

    elseif X(ptr) == X(ptr-1)
        count=1; %計算有幾個重複

        %計算相同類別重複3個以上的情況
        if X(ptr+1) == X(tmp_head)
            %避免1321111這樣的情況出現問題
            if ptr-tmp_head >= 2
                if X(ptr+1) == X(tmp_head+1)
                    while ptr+count<=n_f && X(ptr+count) == X(ptr-1)
                        count = count+1;
                    end
                end
            elseif ptr-tmp_head == 1
                while ptr+count<=n_f && X(ptr+count) == X(ptr-1)
                    count = count+1;
                end
            end
        end

        if count>=2 %相同類別重複3個以上的情況
            ptr = tmp_head + count;
            if ptr < n_f
                Res(partition, 1:count) = X(tmp_head:ptr-1);
            elseif ptr == n_f
                Res(partition, 1:(count+1)) = X(tmp_head:ptr);
            end
        elseif count==1 %沒有重複相同類別的情況
            Xlen = ptr-tmp_head;
            if ptr==length(X)
                Res(partition , 1:Xlen) = X(tmp_head:ptr);
            elseif ptr~=length(X)
                Res(partition , 1:Xlen) = X(tmp_head:(ptr-1));
            end
        end

        tmp_head = ptr; %把暫時的頭移到未計算序列的第一個
        if ptr~=length(X)
            partition = partition + 1;
        end
    end

    if ptr < length(X)
        ptr = ptr+1;
    elseif ptr == length(X)
        break;
    end
end

%% 開始執行演算法

E_T_ptr = 1; %每一組的第一個位置
total_part = total_part + partition;

for round = 1:partition
    %% 計算每一組的長度，記錄在len %%
    Res_R = find(Res(round,:)); 
    len = length(Res_R);
    
    %%
    %判斷此組到底是不是連續序列
    %若=0，表示此組不為連續序列 ； 若=1，表示此組為連續序列 
    conti_flag = 1;
    if len>=2
        for conti_count=2:len
            if Res(round,conti_count) ~= Res(round,conti_count-1)
                conti_flag = 0;
                break;
            end
        end
    end
    
    %% 執行演算法 %%
    % 若不為連續序列 (也就是conti_flag == 0時)、或者並非只有一架飛機的序列，需要執行演算法
    
    %% 讓這組有值的地方變成一
    if conti_flag==0 || len~=1
        H_ = Res(round , 1:len);
        M_ = Res(round , 1:len);
        S_ = Res(round , 1:len);
        for j=1:len
            if( Res(round,j) == 2)
                H_(j) = 1;
                M_(j) = 0;
                S_(j) = 0;
            elseif( Res(round,j) == 3 )
                H_(j) = 0;
                M_(j) = 1;
                S_(j) = 0;
            elseif Res(round,j) == 1 
                H_(j) = 0;
                M_(j) = 0;
                S_(j) = 1;
            end
        end
        
    %%
        E_T_tmp=[]; % 暫時紀錄每一組的E_T
        E_T_tmp(1:len) = E_T( E_T_ptr : E_T_ptr+len-1 ); % E_T_ptr : 每一組的第一台是幾號飛機

        EST_tmp=[]; % 暫時紀錄每一組的EST
        for i=1:len
            EST_tmp(i) = E_T_tmp(i) - E_T_tmp(1);
        end
        EST = EST_tmp';
        %和第一架相減之後變成和第一架飛機相差的進場時間稱為EST

        E_T_ptr = E_T_ptr+len;

        %HG的89.63和 MG的48.68還有HA跟MA的51.75是TA和TB的值
        SG=120.44*S_;
        HG=67.44*H_;
        MG=50.84*M_;
        G=SG+HG+MG;
        
        SA=52.15*S_;
        HA=52.15*H_;
        MA=64.56*M_;
        A=SA+HA+MA;


        t_a=10; %t_a是用來設定飛機可提早和延後的時間(分鐘)
        R=EST-t_a*60; %可提早時間(秒)
        P=EST+t_a*60; %可延後時間(秒)

        t=70;%t是時間窗大小(分鐘)
        range_d=1-(t/80); %80=70+10
        range_i=(t/60)-1; %60=70-10
        v_low=v-v*range_d;
        v_upp=v+v*range_i;
        b_vl=v./(v_low);
        b_vu=v./(v_upp);

        Z1=zeros(1,len);
        Z2=zeros(1,len*len);
        Z3=zeros(1,len-1);
        Z4=zeros(1,len-2);

        I1=ones(1,len);
        I2=ones(1,len-1);

        e1=eye(1,len); %eye(n)特性矩陣(identity matrix)，產生一個nxn單位方矩陣
        e2=zeros(1,len);
        e2(2)=1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %目標函數
        os=0.5;
        f1=os*[kron(EST',e1) I2 Z3 Z3];
        %kron(A,B)函數用於求兩個矩陣的Kronecker積
        %就是矩陣中的每個元素都乘以矩陣B
        f2=(1-os)*t*[Z2 Z3 Z3 I2];
        f=-f1-f2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y1=[];
        y2=[];
        y3=[];
        y4=[];
        y5=[];
        y6=[];
        y7=[];
        y8=[];
        y9=[];
        y10=[];
        Y1=[];
        Y2=[];
        Y3=[];
        Y4=[];
        Y5=[];
        Y6=[];
        Y7=[];
        Y8=[];
        Y9=[];
        Y10=[];

        %A1飛機提早時間限制
        %A2飛機延遲時間限制
        %A3飛機隔離限制
        ee1=zeros(1,len-1);
        for i=1:len-1
            ee1(i)=1;
            ee2=zeros(1,len);
            ee2(i+1)=1;
            ee3=zeros(1,len-1);
            ee3(i)=1;
            ee4=zeros(1,len);
            ee4(i)=1;

            A1=[kron(R',ee2)-kron(EST',e1) -ee1 Z3 Z3];
            A2=[kron(EST',e1)-kron(P',ee2) ee1 Z3 Z3];
            A3=[kron(G,ee4)+kron(A,ee2) -ee3 Z3 Z3];
            size_a=size(A1);
            Y1(i,:)=A1;
            Y2(i,:)=A2;
            Y3(i,:)=A3;
            b1=0;
            b2=0;
            b3=0;
            size_b=size(b1);
            y1(i,:)=b1;
            y2(i,:)=b2;
            y3(i,:)=b3;
        end
        %A4經過time window調整仍須滿足隔離需求
        for i=1:len-2
            eee2=zeros(1,len);
            eee2(i+2)=1;
            eee3=zeros(1,len);
            eee3(i+1)=1;
            eee4=zeros(1,len-1);
            eee4(i)=1;
            eee5=zeros(1,len-1);
            eee5(i+1)=1;

            A4=[-kron(EST',eee2-eee3) eee4 -t*60*(eee5-eee4) Z3];
            size_e=size(A4);
            Y4(i,:)=A4;

            b4=0;
            size_f=size(b4);
            y4(i,:)=b4;
        end
        %A5飛機速度限制(提前)
        %A6飛機速度限制(延後)
        for i=1:len-1
            eee6=zeros(1,len-1);
            eee6(i)=1;

            A5=[Z2 Z3 eee6 Z3];
            A6=[Z2 Z3 -eee6 Z3];
            size_g=size(A5);
            Y5(i,:)=A5;
            Y6(i,:)=A6;

            b5=b_vl(i+1);
            b6=-b_vu(i+1);

            size_h=size(b5);
            y5(i,:)=b5';
            y6(i,:)=b6';
        end


        %A7邊界條件
        A7=[-kron(EST',e2-e1) 1 Z4 -t*60 Z4 Z3];
        size_i=size(A7); %d= size(X);%返回矩陣的行數和列數，保存在d中
        %ex:x=[1,2,3;4,5,6]是一個2*3的矩陣 %d=size(x) d=23
        Y7(1,:)=A7;

        b7=-t*60;
        size_j=size(b7);
        y7(1,:)=b7;

        %提前延後的cost,c1延後c2提前
        c1=50;
        c2=5;   
        for i=1:len-1
            eee7=zeros(1,len-1);
            eee7(i)=1;

            A8=[Z2 Z3 -c1*t*60*eee7 -eee7];
            A9=[Z2 Z3 c2*t*60*eee7 -eee7];
            size_k=size(A8);
            Y8(i,:)=A8;
            Y9(i,:)=A9;

            b8=-c1*60*t;
            b9=c2*60*t;
            size_m=size(b8);
            y8(i,:)=b8;
            y9(i,:)=b9;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %S順序矩陣的限制(每一行每一列加起來都必須等於一)
        s1=repmat(eye(len),[1 len]);%堆疊矩陣3
        s2=reshape(repmat(eye(len),[len 1]),[len len.^2]);
        s3=eye(1,len*len);
        s12=[s1;s2;s3];
        s4=zeros(2*len+1,3*(len-1));
        SS=[s12 s4];
        size_k=size(SS);
        Y10=SS;

        ss=ones(2*len+1,1);
        size_l=size(ss);
        y10=ss;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=[Y1;Y2;Y3;Y4;Y5;Y6;Y7;Y8;Y9;Y10];
        b=[y1;y2;y3;y4;y5;y6;y7;y8;y9;y10];

        vlb=[zeros(1,(length(EST'))^2) 30*ones(1,length(EST')-1) zeros(1,len-1) zeros(1,len-1)];
        vub=[ones(1,(length(EST'))^2) 1000*ones(1,length(EST')-1) 5*ones(1,len-1) 10000*ones(1,len-1)];
        e=[-ones(1,3*(len-1)) -ones(1,len-2) -ones(1,2*(len-1)) -ones(1,1) -ones(1,2*(len-1)) zeros(1,2*len+1)];

        xint=[];
        for i=1:len^2 ;%len^2
            xint(i)=i;
        end
 
        % USE Gurobi
        clear model;
        model.sense(1:8*len-8)=['<'];%sense=The senses of the linear constraints.
        model.sense(8*len-7:10*len-7)=['='];
        %檔案儲存格式：
        %以附檔名lp存檔，例  model.lp
        model.obj=f;
        model.lb=vlb;
        model.ub=vub;
        model.A=sparse(a);
        model.rhs=b;
        model.modelsense='max';
        model.vtype=char([ones(1,length(xint))*'B' ones(1,length(f)-length(xint))*'C']);
        clear params;
        params.threads=4;
        result=gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            %%strcmp 指令用於比較字串內容的異同
            gurobi_time = gurobi_time + result.runtime; %只有跑gurobi的時間
        else %無解的情況，直接跳出for迴圈
            break;
        end
        
        temp_type = []; % 記錄每一組的進場類別結果
        temp_num=[]; % 紀錄各組中每一台最後的進場順序
        
        temp_num = find(result.x); % 暫存result.x第幾個值不等於0
        for i = 1:len
            temp_num(i) = mod(temp_num(i), len);
            if temp_num(i) == 0
                temp_num(i) = len; % 若餘數=0，代表順序沒有變動
            end
            temp_type(temp_num(i)) = Res(round, i);
        end
        res_type = [res_type temp_type];
        temp_num = temp_num + num_tmp;
        res_num_prev = [res_num_prev temp_num(1:len)'];
    else % 是連續序列或者單一序列時，直接在res_order後面接上去即可
        res_type = [res_type Res(round, 1: len)];
        tmp_count=[];
        for i=1:len
            tmp_count(i) = i;
        end
        tmp_count = tmp_count + num_tmp;
        res_num_prev = [res_num_prev tmp_count];
    end
    num_tmp = num_tmp + len;
end
total_time = etime(clock, t0);
avg_time = total_time/total_part;
sprintf('Average computation time for %d aircraft is %f seconds',n_f,avg_time)

for i=1:n_f
   res_num(res_num_prev(i)) = i; 
end

X
res_type
res_num

%% 計算最後降落時間
landing_time = 0;
for i = 2 : n_f
    if res_type(i-1) == 1
        landing_time = landing_time + 120.44;
    elseif res_type(i-1) == 2
        landing_time = landing_time + 67.44;
    else
        landing_time = landing_time + 50.84;
    end
    
    if res_type(i) == 1
        landing_time = landing_time + 52.15;
    elseif res_type(i) == 2
        landing_time = landing_time + 52.15;
    else
        landing_time = landing_time + 64.56;
    end
end
landing_time