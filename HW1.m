close all
clearvars
clc

T_w = 0.01;                                                                 % duration of each window in seconds

[ x , f_s ] = audioread( 'speech.wav' );                                    % loads the speech signal
[ y , ~ ] = audioread( 'piano.wav' );                                       % loads the musical signal

M = T_w * f_s;                                                              % duration of each window in samples

if ( mod( M , 2 ) == 1 )
    M = M + 1;                                                              % if the duration of the window in samples is odd, it makes it even
end

win = hanning( M , 'periodic' );                                            % creates a hanning window

if ( iscola( win , M / 2 , 'ola' ) )
    fprintf( "The window respects the COLA condition\n" );                  % tells us if the COLA condition is satisfied
else
    fprintf( "The window doesn't respect the COLA condition\n" );           % warns us that the COLA condition is not satisfied
end

N_zp = M - mod( length( x ) , M / 2 );                                      % number of zeros to add to the tail of the signal

x = padarray( x , N_zp , 0 , 'post' );                                      % zero-pad the array to match an integer number of windows
y = padarray( y , N_zp , 0 , 'post' );                                      % zero-pad the array to match an integer number of windows

p_x = 46;                                                                   % defines the order of the LPC (number of filter taps) for speech signal
p_y = 46;                                                                   % defines the order of the LPC (number of filter taps) for the musical signal

u = zeros( M , 1 );                                                         % creates an impulse signal to obtain filter impulse responses
u( 1 ) = 1; 

s = zeros( length( x ) , 1 );

for k = 1 : 2 * length( x ) / M - 1
    l_idx = 1 + ( k - 1 ) * M / 2;                                          % set the left index for windowing
    r_idx = M + ( k - 1 ) * M / 2;                                          % set the right index for windowing

    x_win = x( l_idx : r_idx ) .* win;                                      % select the k-th frame of the speech signal (windowing)
    a_x_win = wienerHopf( x_win , p_x );                                    % computes the speech filter coefficients for the k-th frame using Wiener-Hopf equations
    A_x_win = [ 1 , -a_x_win' ];                                            % defines the speech whitening filter for the k-th frame
    h_x_win = filter( 1 , A_x_win , u );                                    % computes the impulse response of the speech shaping filter for the k-th frame
    H_x_win = fft( h_x_win );                                               % computes the frequency response of the speech shaping filter for the k-th frame
    
    y_win = y( l_idx : r_idx ) .* win;                                      % select the k-th frame of the musical signal (windowing)
    a_y_win = wienerHopf( y_win , p_y );                                    % computes the musical filter coefficients for the k-th frame using Wiener-Hopf equations
    A_y_win = [ 1 , -a_y_win' ];                                            % defines the musical whitening filter for the k-th frame
    e_y_win = filter( A_y_win , 1 , y_win );                                % computes the prediction error of the musical signal in the k-th frame
    E_y_win = fft( e_y_win );                                               % computes the spectrum of the prediction error for the k-th frame

    S_win = H_x_win .* E_y_win;                                             % the prediction error of the musical signal is filtered (in frequency domain) with the speech shaping filter in the k-th frame
    s_win = ifft( S_win ) / M;                                              % the k-th frame of the resulting signal is obtained through iFFT

    s( l_idx : r_idx ) = s( l_idx : r_idx ) + s_win;                        % overlap and add process to obtain the entire processed signal
end

figure( 1 );
plot( abs( s ) );



function a = wienerHopf( x , p )

[ r , ~ ] = xcorr( x , x ); % compute the autocorrelation vector of the signal

R_p = toeplitz( r( 1 : p ) ); % creates the autocorrelation matrix of order p

r_p = r( 2 : p + 1 ); % creates the autocorrelation vector of order p

a = R_p \ r_p; % computes the filter coefficients

end
