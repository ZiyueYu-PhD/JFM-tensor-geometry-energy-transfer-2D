function [x_grid,y_grid,eigvec_x,eigvec_y,eig_val,eigvec_angle]=...
    Cal_eigvector_eigvalue_strainandstress(px,py,pvxs,pvys,type,L,d0,filter)


Wn=2*d0/L;
if filter
    pvx=sharpfilt2(pvxs,Wn,0.5,1,0);
    pvy=sharpfilt2(pvys,Wn,0.5,1,0);
    pvxpvx=sharpfilt2(pvxs.*pvxs,Wn,0.5,1,0);
    pvypvy=sharpfilt2(pvys.*pvys,Wn,0.5,1,0);
    pvxpvy=sharpfilt2(pvxs.*pvys,Wn,0.5,1,0);
else
    pvx=pvxs;
    pvy=pvys;
end


[U2x,U2y] = gradient(pvx(:,:),px(1,:),py(:,1));
[V2x,V2y] = gradient(pvy(:,:),px(1,:),py(:,1));

eigvec_angle=zeros(size(U2x));
eigvec_x = zeros(size(U2x));
eigvec_y = zeros(size(U2x));
eig_val = zeros(size(U2x));

% subgrid stress
subt11=sharpfilt2(((pvxs-pvx).^2),Wn,0.5,1,0);
subt12=sharpfilt2(((pvxs-pvx).*(pvys-pvy)),Wn,0.5,1,0);
subt21=sharpfilt2(((pvxs-pvx).*(pvys-pvy)),Wn,0.5,1,0);
subt22=sharpfilt2(((pvys-pvy).^2),Wn,0.5,1,0);

for k=1:size(U2x,1)
    for l=1:size(U2x,2)
        if type==1
            % type=1 means strain
            strain=[U2x(k,l), 0.5*(U2y(k,l)+V2x(k,l));...
                0.5*(U2y(k,l)+V2x(k,l)), V2y(k,l)];
            %         eig_strain(k,l)=max(abs(eig(strain)));

        else
            % type=0 means stress
            % please notice here we use subgrid stress
            % strain=[pvxpvx(k,l)-pvx(k,l)*pvx(k,l), pvxpvy(k,l)-pvx(k,l)*pvy(k,l);...
            %     pvxpvy(k,l)-pvx(k,l)*pvy(k,l), pvypvy(k,l)-pvy(k,l)*pvy(k,l)];
            % strain=strain-1/2*(strain(1,1)+strain(2,2))*eye(2,2);

            stress=[subt11(k,l),subt12(k,l);...
                subt21(k,l),subt22(k,l)];
            strain=stress;
        end

        [V,D] = eig(strain);
        eig_val(k,l) = abs(max(diag(D)));
        eig_val(k,l) = max(diag(D));
        
        % if eig_val(k,l)<0
        %     disp('Eig_val <0')
        %     disp(eig_val(k,l))
        %     disp([k,l])
        %     disp(type)
        %     pause
        % end

        eig_vec = V(:,diag(D)==max(diag(D))); % eigen vector
        % adjust the vector direction
        eigvec_x(k,l) = eig_vec(1)/norm(eig_vec)*sign(eig_vec(1));
        eigvec_y(k,l) = eig_vec(2)/norm(eig_vec)*sign(eig_vec(1));
        % keep original direction
        %         eigvec_x(k,l) = eig_vec(1)/norm(eig_vec);
        %         eigvec_y(k,l) = eig_vec(2)/norm(eig_vec);
        eigvec_angle(k,l) = atand(eigvec_y(k,l)/eigvec_x(k,l));
    end
end
% %%
% figure(2)
x_grid=px; y_grid=py;
% imagesc(x_grid(1,:),y_grid(:,1),eigvec_angle);
% % imagesc(1:1600,1:1600,eigvec_angle);
%
% hold on
% quiver(x_grid,y_grid,eigvec_x,eigvec_y,'color','k');
% pbaspect([1 1 1]);
% set(gca,'Ydir','reverse');
% colorbar;
% title('eigenvector in degree');
%
% axis([0 1024 0 1024])
% set(gca,'YDir','reverse');
% xlb = xlabel('x');
% ylb = ylabel('y');
% set([xlb,ylb],'interpreter','Latex','fontsize',11);
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gcf,'units','pixels','innerposition',[200,200,800,800]);
% set(gca,'DataAspectRatio',[1,1,1]);
% set(gca,'looseInset',[0.015 0.01 0.03 0.01]);
%%
% figure(2)
% imagesc(x_grid,y_grid,eig_val);
% % hold on
% % quiver(x_grid,y_grid,eigvec_x,eigvec_y,'color','k');
% pbaspect([1 1 1]);
% set(gca,'Ydir','reverse');
% colorbar;
% clim([0,0.3])
% title('eigenvalue');
%
% axis([50 950 50 950])
% set(gca,'YDir','reverse');
% xlb = xlabel('x');
% ylb = ylabel('y');
% set([xlb,ylb],'interpreter','Latex','fontsize',11);
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gcf,'units','pixels','innerposition',[200,200,800,800]);
% set(gca,'DataAspectRatio',[1,1,1]);
% set(gca,'looseInset',[0.015 0.01 0.03 0.01]);


% %%
% clear;
% % -=- load files -=-
% mainfolder = '<DATA_ROOT>\set4\';
% filename = 'UV2D_2.mat';
% load([mainfolder,filename],'u','v','u_fluc','v_fluc','mean_u','mean_v','x','y');
% %%
% U = mean_u;
% V = mean_v;
% uu=mean(u_fluc.^2,3);
% vv=mean(v_fluc.^2,3);
% uv=mean(u_fluc.*v_fluc,3); % Raynolds stress
% vu=mean(v_fluc.*u_fluc,3);
% [U2x,U2y]=gradient(U,x,y);
% [V2x,V2y]=gradient(V,x,y);
% % S_ij = 0.5*(U2y+V2x); % off-diagnol entry of the rate of strain of mean flow
%
% % -=- calculate production -=-
% P = -uu.*U2x -uv.*U2y -vu.*V2x - vv.*V2y;
% %%
% imagesc(x,y,P);
% colorbar;
% hold on
% quiver(x_grid,y_grid,eigvec_x,eigvec_y,'color','k');
% set(gca,'Ydir','reverse');
% % colorbar;
% % title('eigenvector in degree');
% %%
% % im = imread([mainfolder,'set8_archive/Basler_acA2040-90um__23513671__20210830_213858735_610.tiff']);
% % imshow(im);
% % hold on
% figure();
% for i = 40
%     quiver(x_grid,y_grid,mean_u,mean_v);
%     hold on
%     quiver(x_grid,y_grid,u_fluc(:,:,i),v_fluc(:,:,i));
%     set(gca,'Ydir','reverse');
%     pause(0.2);
%     hold off
% end
% %%
% u_select = u_fluc(:,:,:);
% v_select = v_fluc(:,:,:);
%
% uu = mean(u_select.^2,3);
% vv = mean(v_select.^2,3);
% uv = mean(u_select.*v_select,3);
% vu = mean(v_select.*u_select,3);
%
% Trace_mean = 0.5*(uu+vv);
%
%
% eigvec_angle=zeros(size(uu));
% for ii=1:size(uu,1)
%     for jj=1:size(uu,2)
%         stress=[uu(ii,jj)-Trace_mean(ii,jj), uv(ii,jj);...
%             vu(ii,jj), vv(ii,jj)-Trace_mean(ii,jj)]; % traceless
%
%         [V,D] = eig(stress);
%         eig_val(ii,jj) = max(diag(D));
%         eig_vec = V(:,diag(D)==max(diag(D)));
%         eigvec_x(ii,jj) = eig_vec(1)/norm(eig_vec)*sign(eig_vec(1));
%         eigvec_y(ii,jj) = eig_vec(2)/norm(eig_vec)*sign(eig_vec(1));
%         eigvec_angle(ii,jj) = atand(eig_vec(2)/eig_vec(1));
%     end
% end
% %%
% figure(1);
% imagesc(x_grid,y_grid,V2x);
% hold on
% quiver(x_grid,y_grid,u_select,v_select,'color','k');
% set(gca,'Ydir','normal');
% colorbar;
% title('dvdx');
% %%
% figure(2)
% imagesc(x_grid,y_grid,eig_val);
% hold on
% quiver(x_grid,y_grid,eigvec_x.*eig_val/max(eig_val(:)),eigvec_y.*eig_val/max(eig_val(:)),'color','k');
% set(gca,'Ydir','reverse');
% colorbar;
% title('eigenvalue and eigenvector');
% %%
% im = imread([mainfolder,'set8_archive/Basler_acA2040-90um__23513671__20210830_213858735_680.tiff']);
% imshow(im);
% hold on
% quiver(x_grid,y_grid,eigvec_x.*eig_val/max(eig_val(:)),eigvec_y.*eig_val/max(eig_val(:)),'color','r');
% % quiver(x_grid,y_grid,u_select,v_select,'color','y');
% set(gca,'Ydir','reverse');
% title('veloctity');
% %%
% im = imread([mainfolder,'set8_archive/Basler_acA2040-90um__23513671__20210830_213858735_90.tiff']);
% imshow(im);
% hold on
% quiver(x_grid,y_grid,eigvec_x,eigvec_y,'color','r');
% % quiver(x_grid,y_grid,u_select,v_select,'color','y');
% set(gca,'Ydir','reverse');
% title('veloctity');

end