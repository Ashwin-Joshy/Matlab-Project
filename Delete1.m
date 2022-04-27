D = uint8('142')
D = de2bi(D,8);
D = reshape(D,[1, size(D,1)*8]);
Data_Key=uint8('111')
E=de2bi(Data_Key,8);
E = reshape(E,[1, size(E,1)*8]);
num_D = length(D); % Find the length of the embedded data D
Encrypt_D = D; % build a container that stores encrypted secrets
%% XOR encryption of the original secret information D according to E
for i=1:num_D
    Encrypt_D(i) = bitxor(D(i),E(i));
end
disp('you are gonna pay:(Encrypt_D)')
disp( Encrypt_D )
K=Encrypt_D
Decrypt_D = K;
num_K=length(K)
for i=1:num_K
    Decrypt_D(i) = bitxor(K(i),E(i)); 
end
isequal(D,Decrypt_D)
%secret1 = reshape(Decrypt_D,[32 8]);
%secret2 = bi2de(secret1);
%secret3 = char(secret2)';
%disp(secret3)