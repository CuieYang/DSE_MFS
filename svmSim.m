function [Yd, Ydreal] = svmSim(svm,Xt)

% ------------------------------------------------------------%

cathe = 10e+6;                 % kernel�����Ԫ�ظ���������
nx = size(svm.x,2);            % ѵ��������
nt = size(Xt,2);               % ����������
block = ceil(nx*nt/cathe);     % �ֿ鴦��
num = ceil(nt/block);          % ÿ�����������

for i = 1:block
    if (i==block)
        index = [(i-1)*num+1:nt];
    else
        index = (i-1)*num+[1:num];
    end

    [Yd(index),Ydreal(index)] = svmSim_block(svm,Xt(:,index));           % �������
end

% ------------------------------------------------------------%

function [Yd, Ydreal] = svmSim_block(svm,Xt);

% �������:
% svm  ֧��������(�ṹ�����)
% the following fields:
%   type - ֧������������  {'svc_c','svc_nu','svm_one_class','svr_epsilon','svr_nu'}
%   ker - �˲���
%       type   - linear :  k(x,y) = x'*y
%                poly   :  k(x,y) = (x'*y+c)^d
%                gauss  :  k(x,y) = exp(-0.5*(norm(x-y)/s)^2)
%                tanh   :  k(x,y) = tanh(g*x'*y+c)
%       degree - Degree d of polynomial kernel (positive scalar).
%       offset - Offset c of polynomial and tanh kernel (scalar, negative for tanh).
%       width  - Width s of Gauss kernel (positive scalar).
%       gamma  - Slope g of the tanh kernel (positive scalar).
%   x - ѵ������
%   y - ѵ��Ŀ��;
%   a - �������ճ���
%
% Xt  ��������,d��n�ľ���,nΪ��������,dΪ����ά��

% �������:
% Yd  �������,1��n�ľ���,nΪ��������,ֵΪ+1��-1

% ------------------------------------------------------------%

type = svm.type;
ker = svm.ker;
X = svm.x;
Y = svm.y;
a = svm.a;

% ------------------------------------------------------------%
% �������

epsilon = 1e-8;                  % ���С�ڴ�ֵ����Ϊ��0
i_sv = find(abs(a)>epsilon);          % ֧�������±�

switch type
    case 'svc_c',
        
        tmp = (a.*Y)*kernel(ker,X,X(:,i_sv));          % ������
        b = Y(i_sv)-tmp;
        b = mean(b);
        tmp =  (a.*Y)*kernel(ker,X,Xt);
        tmp = tmp+b;
        Yd = sign(tmp);
        
    case 'svc_nu', 
        %------------------------------------%
        % �� 'svc_c' �����ͬ

        tmp = (a.*Y)*kernel(ker,X,X(:,i_sv));          % ������
        b = Y(i_sv)-tmp;
        b = mean(b);
        tmp =  (a.*Y)*kernel(ker,X,Xt);
        Yd = sign(tmp+b);
        
    case 'svm_one_class',
        
        n_sv = length(i_sv);
        tmp1 = zeros(n_sv,1);
        for i = 1:n_sv
            index = i_sv(i);
            tmp1(i) = kernel(ker,X(:,index),X(:,index));
        end

        tmp2 = 2*a*kernel(ker,X,X(:,i_sv));           % ������
        tmp3 = sum(sum(a'*a.*kernel(ker,X,X)));    

        R_square = tmp1-tmp2'+tmp3;
        R_square = mean(R_square);                       % ����뾶 R^2 (������֧��������ƽ���Ľ��)

        nt = size(Xt,2);                  % ����������

        tmp4 = zeros(nt,1);               % ������
        for i = 1:nt
            tmp4(i) = kernel(ker,Xt(:,i),Xt(:,i));
        end
    
        tmp5 = 2*a*kernel(ker,X,Xt);                % ������
        Ydreal = tmp4-tmp5'+tmp3-R_square;
        Yd = sign(Ydreal);
        

    case 'svr_epsilon',
        
        tmp = a*kernel(ker,X,X(:,i_sv));   % ������
        b = Y(i_sv)-tmp;                    % ���Ų�һ��,���ߺ����Ͳ�һ��,ʵ������һ����!
        %b = Y(i_sv)+tmp;
        b = mean(b);

        tmp =  a*kernel(ker,X,Xt);         % ���Ų�һ��,���ߺ����Ͳ�һ��,ʵ������һ����!
        %tmp =  -a*kernel(ker,X,Xt);
        Yd = (tmp+b);        
        
    case 'svr_nu',
        %------------------------------------%
        % ��'svr_epsilon' �����ͬ
        
        tmp = a*kernel(ker,X,X(:,i_sv));   % ������
        b = Y(i_sv)-tmp;                    % ���Ų�һ��,���ߺ����Ͳ�һ��,ʵ������һ����!
        %b = Y(i_sv)+tmp;
        b = mean(b);

        tmp =  a*kernel(ker,X,Xt);         % ���Ų�һ��,���ߺ����Ͳ�һ��,ʵ������һ����!
        %tmp =  -a*kernel(ker,X,Xt);
        Yd = (tmp+b);        
        
    otherwise,
end