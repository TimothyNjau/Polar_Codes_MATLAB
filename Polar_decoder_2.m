R = [1 2 3 5 9 4 6 10 7 11 13 8 12 14 15 16]; %reliability sequence of N=16.Ordered from worst to best.

N = 16; %N bits
n = log2(N);
K = 8; %K bits
R1 = R(R<=N); % the reliability sequence of N=8
EbNodb = 4;
Rate = K/N; %transmission rate
EbNo = 10^(EbNodb/10);
sigma = sqrt(1/(2*Rate*EbNo));

%polar encoder
F = R1(1:N-K);% frozen positions in sequence: R(1:N-K). These are the initial positions
% of the reliability channel.
% messsage positions in sequence : R(N-K+1:end)

messg = ([1 1 0 1 0 1 0 1]);%randi([0 1],1,K); % generates random message bits of K size

y = zeros(1,N); %assign zero to frozen positions.
y(R1(N-K+1:end)) = messg; %assign message bits to reliable positions.

m = 1; %number of bits combined

for d = n-1:-1:0
    for i = 1:2*m:N
        a = y(i:i+m-1); %first part
        b = y(i+m:i+2*m-1); %second part
        y(i:i+2*m-1) = [mod(a+b,2) b]; %combining
    end
    m = m * 2;
end

%transmission of polar bits
wordcod = y; %final transmitted coded bits
s = 1-2*wordcod; %BPSK bit to symbol conversion
r = s + sigma * randn(1,N); %receiver channel with added AWGN


%successive cancellation decoder
node = 0; depth = 0; %start of root
L = zeros(n+1, N); %beliefs
ns = zeros(1,2*N-1); %node state vector
L(1,:) = r; %belief of root

f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).* min(abs(a),abs(b)); %minsum
g = @(a,b,c) b+(1-2*c).*a; %g function

done = 0; %decoder is finished or not
while (done==0) %traverse till all bits are decoded
    %determine if node is leaf or not
    if depth == n
        if any (F==(node+1)) %if node is frozen
            ucap(n+1,node+1) = 0;
        else
            ucap(n+1,node+1) = 1;
        end
        node = floor(node/2); depth = depth-1;
    else
        %for nonleaf
        npos = (2^depth - 1) + node + 1;%position of node in state vector
        if ns(npos) == 0
            temp = 2^(n-depth); %temp storage
            Ln = L(depth+1, temp*node+1:temp*(node+1)); %incoming beliefs
            
            %split beliefs into two
            a = Ln(1:temp/2);
            b = Ln((temp/2) +1:end);
            
            node = node*2; depth= depth + 1; %next node: left child
            temp = temp/2; %incoming belief length for left child
            
            L(depth+1, temp*node+1:temp*(node+1)) = f(a,b); %minsum and storage
            ns(npos) = 1;
            
        else
            if ns(npos) == 1
                temp = 2^(n-depth); %temp storage
                Ln = L(depth+1, temp*node+1:temp*(node+1)); %incoming beliefs
                
                %split beliefs into two
                a = Ln(1:temp/2);
                b = Ln((temp/2) +1:end);
                
                lnode = 2*node; ldepth = depth + 1;
                ltemp = 2^(n-ldepth);
                ucapn = ucap(ldepth + 1, ltemp*lnode+1:ltemp*(lnode+1)); %incoming declarations from left child
                node = node*2 + 1; depth= depth + 1; %next node: right child
                temp = temp/2; %incoming belief length for right child
                
                L(depth+1, temp*node+1:temp*(node+1)) = g(a,b,ucapn); %g and storage
                ns(npos) = 2;
            else
                temp = 2^(n-depth);
                lnode = 2*node; rnode = 2*node; cdepth = depth + 1;
                ctemp = temp/2;
                ucapl = ucap(cdepth + 1, ctemp*lnode+1:ctemp*(lnode+1));
                ucapr = ucap(cdepth + 1, ctemp*lnode+1:ctemp*(rnode+1));
                ucap(depth+1, temp*node+1:temp*(node+1)) = [mod(ucapl + ucapr,2) ucapr];
                node = floor(node/2); depth = depth -1;
                
            end
        end
    end
    
end










