%單跑道模擬

clc
clear all

format short
op_sch_all=[];
op_ori_all=[];
Time_vartotal_all=[];
Time_no_all=[];
Time_final_all=[];
FCvar_all=[];

for nua=1:1    % 隨機?筆
    ptw=60    %  width of controlled time windows
    
    %% 初始值
    n_f=10;    % Number of aircraft
%     X= randint(1,n_f,[1 3]) ;
    X = [1,3,2,3,1,2,3,1,1,1]; %別動他 亂動揍你
    % H_=floor(rand(1,n_f)*2)
    % M_=abs(floor(H_-0.6))
    %E_T = sort(ceil((n_f*105)*rand(1,n_f))); %隨機到場時間
    E_T = [55,114,429,474,483,540,579,736,846,916]
    %E_T = [18 91 127 202 246 453 487 594 598 601 774 987 1084 1248 1427] ;
    EST=[E_T-E_T(1)]';    % 設定 proceding aircraft 到場時間為零
    % v=[ceil(220*rand(1,n_f))+369]  % 速度隨機 %巡航速度
    %v=[507 485 488 512 523 529 509 530 508 506 497 502 505 484 525];
    v=[418,471,543,469,383,443,580,532,481,413];
    lest = length(EST) ;
    partition = 1;
    par_ptr = 1;
    
    
    %% 設定Map
    instr = '';
    %     input = [1 1 2 2 1 3 4 2 1 2  3  4  1  2  3  1  1  3  3  2];
    %            1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
    for i=1:length(X)
        instr = strcat(instr, int2str(X(i)));
    end
    
    keySet = [237 158 264 290 291 317 318 343 344 370 371 396 397 423 449];
    
    cont237 = [1111 2111 3111 3211 3311 3321 3331 3332 3333]; %1
    cont158 = [111 211 311 321 331 332 333]; %2
    cont264 = [2211 3221 3322]; %2
    cont290 = [3231 3232 3233 3323 2332 2333 2321 2331 2311]; %3
    cont291 = [3222 2221]; %4
    cont317 = [3212 3223 3312 3121 3112 2322 2233 2231 2232 2112 2121 1211 1112 1121]; %5
    cont318 = [2222]; %6
    cont343 = [3133 3313 3132 3131 3213 3113 2323 2131 2132 2133 2113 1331 1332 1333 1311 1321 1133 1131 1132 1113]; %7
    cont344 = [3122 2223 2212 2122 1221 1122]; %8
    cont370 = [3123 2312 2213 2123 1233 1322 1231 1232 1123]; %9
    cont371 = [1222]; %10
    cont396 = [2313 1323]; %11
    cont397 = [1212 1223]; %12
    cont423 = [1312 1213]; %13
    cont449 = [1313]; %14
    
    valueSet = {cont237, cont158, cont264, cont290, cont291, cont317, cont318, cont343, cont344, cont370, cont371, cont396, cont397, cont423, cont449};
    M = containers.Map(keySet, valueSet, 'UniformValues', false)
    
    same_pos = []; % same_pos拿來記錄符合14種類型的位置, 範例如下
    % 相同位置  第幾類別  花費時間  機型順序
    % 3         4         291       3222        -> 在X的第3~6台飛機是3222, 他是第4個時間類別, 費時291單位時間
    % 7         12        397       1223        -> 在X的第7~10台飛機是1223, 他是第12個時間類別, 費時397單位時間
    % 12        14        449       1313        -> 在X的第12~15台飛機是1313, 他是第14個時間類別, 費時449單位時間
    
    row = 1; % 初始same_pos的列數
    
    for i=1:2 % 因為有14種時間類別, 所以要跑14次for
        strptr = 1; % 當作instr的暫時指標
        repeat_flag = 0; % 拿來判斷費時較長的時間類別是否重疊到費時短的時間類別
        for j=1:length(valueSet{1,i}) % value{1,i}紀錄每一個時間類別的所有飛機順序組合
            tmp_in = instr(strptr:length(instr)); % 同一個時間類別判斷過的部分不用再判斷一次, 所以重設一個暫時的input
            test_same = findstr(tmp_in, int2str(valueSet{1,i}(j))); % test_same紀錄tmp_in在第幾個位置判斷出相同順序
            if isempty(test_same) % 若沒有判斷出有相同順序, 直接開始判斷下一個時間類別
                continue;
            else % 若有判斷出有相同順序
                % 判斷是否有與較短時間花費的類別重疊的部分
                for k=1:row
                    if i==1 % 在判斷第一類別的時候
                        if isempty(same_pos)==0 && test_same >= same_pos(k,1)-4 && test_same <= same_pos(k,1)+3
                            repeat_flag = 1; % 若有重疊的部分, flag設為1
                            break;
                        end
                    else % 在判斷第二類別的時候
                        if isempty(same_pos)==0 && test_same >= same_pos(k,1)-3 && test_same <= same_pos(k,1)+2
                            repeat_flag = 1; % 若有重疊的部分, flag設為1
                            break;
                        end
                    end
                end
                if repeat_flag == 0 % 若flag=0, 代表沒有重複的部分
                    tmp_sort = [test_same i keySet(i) valueSet{1,i}(j)];
                    same_pos = cat(1, same_pos, tmp_sort); % 把新發現的順序以及類別記錄到same_pos
                    if i==1 % 判斷第一類別時
                        strptr = test_same + 4; % 更新strptr的位置
                    else % 判斷第二類別時
                        strptr = test_same + 3; % 更新strptr的位置
                    end
                end
                mat_size = size(same_pos);
                row = mat_size(1); % 更新same_pos的列數
            end
        end
    end
    
    %%
    H_ = X;
    M_ = X;
    S_ = X;
    for i=1:n_f
        if( X(i) == 2)
            H_(i) = 1;
            M_(i) = 0;
            S_(i) = 0;
        elseif( X(i) == 3 )
            H_(i) = 0;
            M_(i) = 1;
            S_(i) = 0;
        elseif( X(i) == 1 )
            H_(i) = 0;
            M_(i) = 0;
            S_(i) = 1;
        end
    end
    
    A_=S_'*S_*79;  %SS
    B_=H_'*S_*159; %SH
    C_=M_'*S_*185; %SM
    D_=S_'*H_*79;  %HS
    E_=H_'*H_*106; %HH
    F_=M_'*H_*132; %HM
    G_=S_'*M_*79;  %MS
    H_=H_'*M_*79;  %MH
    I_=M_'*M_*79;  %MM
    
    
    T_sep=A_+B_+C_+D_+E_+F_+G_+H_+I_
    %The separation time matrix of the relationship between aircraft
    
    
    %%
    T_fake=zeros(lest,lest);
    T_fake2=triu(ones(lest,lest),1);      % 逆向區間差 %wj??
    for num=1:1:length(EST-1)
        T_fake=T_fake+triu(ones(lest,lest),num);        %  逆向改變次數 %zj??
    end
    %triu 上三角矩陣
    
    %T_fake=triu(ones(length(EST),length(EST)),1);
    %T_fake2=zeros(length(EST),length(EST));
    %for num=1:1:length(EST-1)
    %T_fake2=T_fake2+triu(ones(length(EST),length(EST)),num);
    %end
    
    S_star=diag(ones(lest-1,1),1);
    %控制區段時間
    %diag  %X = diag(v,k)
    %以向量v的元素作為矩陣X的第k條對角線元素，當k=0時，v為X的主對角線
    %當k>0時，v為上方第k條對角線；當k<0時，v為下方第k條對角線。
    va=100;
    t_a=4.28;  %可前後調整的時間
    %%
    %R=(EST-(vt/ptg)*1'-T);  % 提前最多
    %P=(EST+(vt/ptg)*1'-T);  % 延後最多
    R=EST-t_a*60;              % 最早到場時間
    P=EST+t_a*60;              % 最晚到場時間
    
    e=[];
    
    I=diag(ones(1,lest),0);    % 單位矩陣
    %diag 對角線矩陣(diagonal matrix)，僅對角線項有值。
    
    I=eye(lest); %eye(n) identity matrix，產生一個nxn單位方矩陣。
    II=eye(lest-1);
    III=eye(lest-1+lest-1+lest-1);
    
    k=[];
    for i=1:1:lest
        k=[k kron(eye(lest),I(i,:)')];          %轉換系數
    end
    
    BB=kron(eye(lest),ones(1,lest));         %  S row e相加等於1
    G=kron(ones(1,lest),eye(lest));          %  S column 相家等於1
    
    gama=0.9;         % weight value
    Delta=0.95;
    %% miqqQG in Matlab 這裡開始是TOMLAB
    % miqqQG is a small example problem for defining and solving
    % mixed-integer quadratic programming problems with quadratic constraints
    % using the TOMLAB format.
    
    Name = 'MIQQ Test Problem 1';
    f_Low = -1E5;
    x_opt = [];
    f_opt = [];
    IntVars = logical([ones(1,lest^2) ones(1,lest-1) ones(1,lest-1) ones(1,lest-1)]); % 3rd variable is integer valued
    
    %logical : Convert numeric values to logicals
    
    
    F=zeros(lest^2+lest-1+lest-1+lest-1,lest^2+lest-1+lest-1+lest-1);            %    OBJECTIVE FUNCTION   (二次式)
    A=[[[BB;G] zeros(lest*2,lest-1+lest-1+lest-1)];[zeros(lest-1+lest-1+lest-1,lest^2) III]];     %限制式合併
    b_L=[ones(lest*2,1);79*ones(lest-1,1);zeros(lest-1,1);zeros(lest-1,1)];             %  變數下限
    b_U=[ones(lest*2,1);1000*ones(lest-1,1);va*ones(lest-1,1);va*ones(lest-1,1)];     %  變數上限
    
    c=Delta*((gama)*([kron(EST',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]'+[zeros(1,lest^2) ones(1,lest-1) zeros(1,lest-1) zeros(1,lest-1)]')+(1-gama)*([zeros(1,lest^2+lest-1) ones(1,lest-1) zeros(1,lest-1)]'))+(1-Delta)*[zeros(1,lest^2+lest-1+lest-1) ones(1,lest-1)]';   %    OBJECTIVE FUNCTION   (一次式)
    
    x_0=zeros(1,(lest)^2+lest-1+lest-1+lest-1)';
    x_L=[zeros(1,(lest)^2) 79*ones(1,lest-1) zeros(1,lest-1) zeros(1,lest-1)]';      %上限
    x_U=[ones(1,(lest)^2) 1000*ones(1,lest-1) va*ones(1,lest-1) va*ones(1,lest-1)]'; % 下限
    x_min=[zeros(1,(lest)^2) 79*ones(1,lest-1) zeros(1,lest-1) zeros(1,lest-1)]';          % 變數最小值
    x_max=[ones(1,(lest)^2) 1000*ones(1,lest-1) va*ones(1,lest-1) va*ones(1,lest-1)]';      % 變數最大值
    
    
    % Adding quadratic constraints
    clear qc
    
    % Q should be spare matrix
    
    on=zeros(1,lest^2+lest-1+lest-1+lest-1);
    
    for j=2:1:lest
        
        on(lest^2+j-1)=1;  % The last constraints would be used  L1 , L1+L2, L1+L2+L3...
        %sum(L)>=e_j*S*R-e_1*S*EST
        Fir_cs=zeros(lest^2+lest-1+lest-1+lest-1,lest^2+lest-1+lest-1+lest-1);
        Fir_cs(1,2)=1;
        Fir_cs(2,1)=-1;
        
        qc(j-1).Q=sparse(Fir_cs);
        qc(j-1).a=[kron(R',I(j,:)) zeros(1,lest-1+lest-1+lest-1)]'-[kron(EST',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]'-on';
        qc(j-1).r_U=0;
        
        %sum(L)<=Q(j)-e_1*S*EST
        Sec_cs=zeros(lest^2+lest-1+lest-1+lest-1,lest^2+lest-1+lest-1+lest-1);
        Sec_cs(1,3)=1;
        Sec_cs(3,1)=-1;
        qc(j-1+lest-1).Q=sparse(Sec_cs);
        qc(j-1+lest-1).a=-[kron(P',I(j,:)) zeros(1,lest-1+lest-1+lest-1)]'+[kron(EST',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]'+on';
        qc(j-1+lest-1).r_U=0;
        
        
        
        % Lj>=e_j*S^{*}*S*S^{T}*e_j
        QCST=[kron(eye(lest),(I(j-1,:)*S_star)');zeros(lest-1+lest-1+lest-1,lest)]*[kron(I(j-1,:),T_sep)*k zeros(lest,lest-1+lest-1+lest-1)];
        qc(j-1+2*(lest-1)).Q=sparse(QCST+QCST');              % it should be sparse matrix
        qc(j-1+2*(lest-1)).a=-2*[zeros(1,lest^2) II(j-1,:) zeros(1,lest-1) zeros(1,lest-1)]';
        qc(j-1+2*(lest-1)).r_U=0;
        
        %Lj>=e_j*S*R-(e1*S*EST+sum(l))
        For_cs=zeros(lest^2+lest-1+lest-1+lest-1,lest^2+lest-1+lest-1+lest-1);
        For_cs(1,2)=1;
        For_cs(2,1)=-1;
        qc(j-1+3*(lest-1)).Q=sparse(For_cs);
        qc(j-1+3*(lest-1)).a=-[kron(EST',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]'+[kron(R',I(j,:)*S_star) zeros(1,lest-1+lest-1+lest-1)]'-[zeros(1,lest^2) II(j-1,:) zeros(1,lest-1) zeros(1,lest-1)]'-on';
        qc(j-1+3*(lest-1)).r_U=0;
        
        
        
        % w>=t_fake   逆向區間差
        QCSF=[kron(eye(lest),(I(j-1,:)*S_star)');zeros(lest-1+lest-1+lest-1,lest)]*[kron(I(j-1,:),T_fake)*k zeros(lest,lest-1+lest-1+lest-1)];
        qc(j-1+4*(lest-1)).Q=sparse(QCSF+QCSF');              % it should be sparse matrix
        qc(j-1+4*(lest-1)).a=-2*[zeros(1,lest^2) zeros(1,lest-1) II(j-1,:) zeros(1,lest-1)]';
        qc(j-1+4*(lest-1)).r_U=0;
        
        %z>t_fake2   逆向改變次數
        QCSFT=[kron(eye(lest),(I(j-1,:)*S_star)');zeros(lest-1+lest-1+lest-1,lest)]*[kron(I(j-1,:),T_fake2)*k zeros(lest,lest-1+lest-1+lest-1)];
        qc(j-1+5*(lest-1)).Q=sparse(QCSFT+QCSFT');
        qc(j-1+5*(lest-1)).a=-2*[zeros(1,lest^2) zeros(1,lest-1) zeros(1,lest-1) II(j-1,:)]';
        qc(j-1+5*(lest-1)).r_U=0;
        
    end
    
    %R1<    < P1  第一架飛機的範圍
    qc(j-1+1+5*(lest-1)).Q=sparse(Fir_cs);
    qc(j-1+1+5*(lest-1)).a=[kron(R',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]'-[kron(EST',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]';
    qc(j-1+1+5*(lest-1)).r_U=0;
    
    qc(j-1+2+5*(lest-1)).Q=sparse(Sec_cs);    %sum
    qc(j-1+2+5*(lest-1)).a=-[kron(P',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]'+[kron(EST',I(1,:)) zeros(1,lest-1+lest-1+lest-1)]';
    qc(j-1+2+5*(lest-1)).r_U=0;
    
    
    
    
    Prob = miqqAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, qc,...
        IntVars, [], [], [],...
        Name, [], [],...
        x_min, x_max, f_opt, x_opt);
    
    Result = tomRun('cplex', Prob, 1)
    
    S=reshape((Result.x_k(1:lest^2)),lest,lest)          %找出S矩陣
    
    Separation_time=Result.x_k(lest^2+1:lest^2+lest-1)      %算出隔離時間
    
    %%
    
    %計算排序後每一架的到場時間(未考慮變動)
    tisu=[];
    for po=1:1:length(Separation_time)
        tisu=[tisu sum(Separation_time(1:po))];
    end
    
    %計算排序後每一家的到場時間(考慮變動)
    tisu_ori=[];
    wm=S*EST;
    tisu=[wm(1) tisu];
    Totalmove_sch=sum(abs((S*EST)'-tisu))
    
    
    
    % FCFS 到場時間
    FCSeparation_time=[];
    for i=1:1:(lest-1)
        FCSeparation_time=[FCSeparation_time I(i,:)*S_star*eye(lest)*T_sep*eye(lest)*I(i,:)']
    end
    
    
    % FCFS 的變動時間
    FCtisu=[];
    for po=1:1:length(FCSeparation_time)
        FCtisu=[FCtisu sum(FCSeparation_time(1:po))];
    end
    FCvar=[EST(1) EST(1)+FCtisu]-EST'      % FCFS 與預寄到場時間 的差距(變動量)
    
    
    
    sep_ori=[];
    for j=1:1:lest-1
        sep_ori=[sep_ori I(j,:)*S_star*eye(lest)*T_sep*eye(lest)'*I(j,:)'];           %  FCFS 相鄰飛機的隔離時間
    end
    
    for j=1:1:length(sep_ori)
        tisu_ori=[tisu_ori sum(sep_ori(1:j))];               %  FCFS 每一架飛機的到場時間
    end
    tisu_ori=[EST(1) tisu_ori];
    
    
    
    
    
    
    eval(['op_sch_',num2str(nua),'=tisu(lest)'])                       % 排序後operation time
    eval(['op_ori_',num2str(nua),'=tisu_ori(lest)'])              % 排序前operation time
    
    
    
    v=(S*v')';
    range_i=(ptw/55.72)-1                % increase                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    range_d=1-(ptw/64.28)                   % decrease                %藉由controlled time windows 的長度來調整 arrival time windows
    v_low=v-v*range_d;                 %最低速度
    v_upp=v+v*range_i;                 %最高速度
    sep_t=Separation_time'/60           % 相鄰所需隔離分鐘數
    
    arr_con_t=[];
    
    %q=S*EST_moti        % two runwat
    q=S*EST                % one runway
    for j=1:1:lest-1
        q_c=q(j+1)-q(j);
        arr_con_t=[arr_con_t q_c/60];              %進入控制區段的時間差距
    end
    
    
    b_t=-(sep_t-arr_con_t)/ptw         % (l_j    -   a_j,j+1    )/ ptw
    
    b_vu=-v./(v_upp);   % 速度上限
    b_vl=v./(v_low);    % 速度下限
    
    bv=[b_vu;b_vl];
    bv=reshape(bv,length(b_vu)+length(b_vl),1)
    B=[b_t';bv];                                 %  x 值的範圍
    
    A1=zeros(lest-1,lest);
    for jj=1:1:lest-1
        A1(jj,jj)=1;
        A1(jj,jj+1)=-1;                %x_j+1 - x_j
    end
    A2=kron(eye(lest),[-1 1]');        %-x_j , x_j
    Aq=[A1;A2];
    
    
    
    
    H=2*eye(length(b_vu));
    f=-2*ones(1,length(b_vu))'
    
    lb = zeros(length(b_vu),1);
    [x,fval,exitflag,output,lambda] = quadprog(H,f,Aq,B,[],[],lb)                        %利用二次規劃法求得x
    
    
    x_n=x';
    
    
    %x_n=x'
    
    
    
    
    %    v_n=[S*reshape(repmat(v,2,1),1,length(EST)*2)']'./x          %two runway
    v_n=v./x_n        %one runway
    
    %  (([S*reshape(repmat(v,2,1),1,length(EST)*2)']'./v_n)*70-70)*60            % Two runway
    %Time_var_1=sum(abs((([S*v']'./v_n)*(70/60)-(70/60))*60*60));      % one runway
    
    
    
    
    %%
    %
    % <<FILENAME.PNG>>
    %
    
    
    op_sch_all=eval(['[op_sch_all op_sch_',num2str(nua),']'])                % 排序後 optimum的所有operation time
    op_ori_all=eval(['[op_ori_all op_ori_',num2str(nua),']'])                % 排序後 FCFS   的所有operation time
    
    
    
    
    
    eval(['FCvar_',num2str(nua),'=sum(abs(FCvar))'])                                                           % 每一架飛機FCFS與EST的差距
    eval(['Time_var_',num2str(nua),'=(([v]./v_n)*(ptw/60)-(ptw/60))*60*60'])                       %排序完移動 (與EST的差距 )
    eval(['Time_vartotal_',num2str(nua),'=sum(abs((([v]./v_n)*(ptw/60)-(ptw/60))*60*60))'])                                   %排序完移動和
    eval(['Time_arrtime_',num2str(nua),'=(([v]./v_n)*(ptw/60)-(ptw/60))*3600+(S*EST)'''])                  %最終到場時間
    
    eval(['Time_no_',num2str(nua),'=sum(abs(tisu-(S*EST)''))'])                                                             %排序完不移動和
    
    
    
    
    Time_vartotal_all=eval(['[Time_vartotal_all Time_vartotal_',num2str(nua),']'])                   %排序完移動(全部飛機)時間變化
    Time_no_all=eval(['[Time_no_all Time_no_',num2str(nua),']'])                                     %排序完不移動(全部飛機)時間變化
    FCvar_all=eval(['[FCvar_all FCvar_',num2str(nua),']'])                                           %FCFS 的變動
    
end

same_pos

%% 紀錄最後序列
Res = [];
Res_type = [];
for i=1:n_f
    tmp = find(S(:,i));
    Res(tmp) = i;
    Res_type(tmp) = X(i);
end

Res
Res_type