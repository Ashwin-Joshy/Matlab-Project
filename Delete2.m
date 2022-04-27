% Dummy secret key:
secret_key = 12345;

% String to encrypt:
D = double('Very secret string asdadaaa');
D = de2bi(D,8).';
D = D(:).';

% Encryption binary array
rand('seed',secret_key); % set the seed
E = round(rand(1,numel(D))*1);

% crypted string
crypted = bitxor(D,E);
disp("encrypted")
disp(crypted)
crypted1=bitxor(crypted,E)
% Decrypted string
rand('seed',secret_key); % set the seed
E = round(rand(1,numel(crypted1))*1);
decrypted = char(sum(reshape(bitxor(crypted,E),8,[]).*(2.^(0:7)).'))
% decrypted = 'Very secret string'
disp("decrypted")
disp(decrypted)