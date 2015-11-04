% % % on 127.0.0.1
% u=udp('127.0.0.1', 9090,'LocalPort', 9091);
% fopen(u)
% fwrite(u, 'This is a test');
% fclose(u)

tic
for i = 1:100
    
    u=udp('127.0.0.1', 9090,'LocalPort', 9091);
    fopen(u)
    fwrite(u, num2str(i));
    fclose(u)
end % END FOR
toc