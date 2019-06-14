function make_movie(handle,rho)

%% make movie
fps = 24;
fname = sprintf('LorenzQpot_rho%.2f.avi',rho);
writerObj = VideoWriter(fname); % Name it.
writerObj.FrameRate = fps; % How many frames per second.
open(writerObj); 

if rho < 1
    avec = [-8,8,-8,8,-4,6];
else
    if rho > 10 & rho < 13
        avec = [-16,16,-16,16,0,20];
    else if rho > 13 & rho < 16
            avec = [-20,20,-20,20,-5,30];
        else if  rho > 24.06 & rho < 24.74
            avec = [-35,35,-35,35,-10,50];
            else if rho == 20
                     avec = [-25,25,-25,25,0,35];
                else
                    avec = [-100,100,-100,100,-5,180];
                end
            end
        end
    end
end
axis(avec);
frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
writeVideo(writerObj, frame);

for a = 1 : 360
    rotate(handle,[0 0 1],1);
    axis(avec);
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end

close(writerObj); % Saves the movie.

end