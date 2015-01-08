function showimage(m, cmap);
% plots an image of matrix M with optional colormap CMAP

% the main idea of this function is that a the image is drawn
% from scatch the first time this function is called only.
% during subsequent function calls, the image is 'refreshed' by
% replacing the image data in the image's own memory.

% We detect first or subsequent calls by finding out whether the variable
% PLOT_HANDLE is defined.
% Because normally local (function scope) variables are not preserved
% throughout subsequent calls, we make it a `peristent' variable.
% When a variable is made peristent, it is created at the same time. In this
% case, the initial value is set to [] (the empty matrix). We can exploit
% this behaviour of matlab for testing if this was the first call to the
% function (PLOT_HANLDE is []) or a subsequent call (PLOT_HANDLE is not [])

if nargin<2,
   cmap = colormap;
end;

persistent PLOT_HANDLE;				% label PLOT_HANDLE as peristent

if isempty(PLOT_HANDLE),				% PLOTHANDLE not defined yet?
   % first call.						%   no, first call
   figure(1);							%   create new figure
   clf;									%   clear it
   set(gcf,'doublebuffer','on');		%   make it fast
   PLOT_HANDLE = image(m);				% 	create a new image object
   if ~isempty(cmap),					%   optional colormap given?
      colormap(cmap);					%     yes, use it.
   end;
   colorbar;
else,									% PLOT_HANDLE does exist
   %   so we have 2nd 3th etc call
   set(PLOT_HANDLE, 'CData', m);		%   referesh image data
   drawnow;								%   and force redrawing NOW
end;


