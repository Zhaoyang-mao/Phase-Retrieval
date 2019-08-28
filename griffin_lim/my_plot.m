function my_plot(x,x_hat_gla)
subplot(2,2,1);
stem(x,'r');
hold on
stem(x_hat_gla,'b');
legend('Source','Recover');
hold off
subplot(2,2,2);
plot(x,'r-');
hold on
plot(x_hat_gla,'b.-');
hold off
legend('Source','Recover');
subplot(2,2,3);
spectrogram(x);
title('Source');
subplot(2,2,4);
spectrogram(x_hat_gla);
title('Recover');
end
