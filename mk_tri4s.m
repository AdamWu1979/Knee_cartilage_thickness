function [tria,xyza,nid,nt] = mk_tri4s(data,kid,idig,iplt);
%MK_TRI4S Makes a triangular mesh by using the ordered slice data from
%         the axially digitized MRI patella bone data for seven specific
%         patellas.
%
%         [TRIA,XYZA] = MK_TRI4S(DAT,KID,IDIG) given a cell array
%         containing three (3) columns matrices with slice coordinate
%         point data, DAT, a three or four character knee identifier,
%         KID, and an integer 1 or 2 for the first (1) or second (2)
%         digitization, returns the three (3) column triangle
%         connectivity matrix, TRIA.  The coordinates for the patella
%         including an additional end point are returned in a three (3)
%         column matrix, XYZA.
%
%         [TRIA,XYZA,NID,NT] = MK_TRI4S(DAT,KID,IDIG) returns the index
%         to the additional end point, NID, and the number of returned
%         triangles, NT.
%
%         [TRIA,XYZA] = MK_TRI4S(DAT,KID,IDIG,IPLT) given a logical true
%         (nonzero), IPLT, will plot the lower patella mesh and the
%         complete mesh in two separate figures.
%
%         NOTES:  1.  For axially digitized bone data for knees:  4_R,
%                 56_R, 57_L and 64_R in the first digitization and 
%                 56_R, 57_L, 64_R and 66_R in the second digitization.
%
%                 2.  Each slice coordinate data matrix must correspond
%                 to one index into the cell array DATA.
%
%                 3.  Two separate functions are used to mesh the
%                 patellas in two parts (upper and lower part).  See
%                 mk_tri4i.m and mk_tri4p2.m for more information on
%                 the methods of triangulations.
%
%                 4.  The M-files in_tri2d.m, isect.m, mk_tri4i.m,
%                 mk_tri4p2.m, nod2tri.m, plane_fit.m, sl_info.m,
%                 tri_fix2.m, tri_norm.m and tri_shr.m must be in the
%                 current path or directory.
%
%         01-Feb-2016 * Mack Gardner-Morse
%
%         10-Feb-2016 * Mack Gardner-Morse * Updated to work for both
%                                            digitization.
%
%         04-Mar-2016 * Mack Gardner-Morse * Updated to work for knee
%                                            4_R.
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)||isempty(iplt)
  iplt = false;
end
%
if (nargin<3)||isempty(idig)
  error(' *** ERROR in MK_TRI4S:  No digitization number!');
end
%
if idig~=1&&idig~=2
  error(' *** ERROR in MK_TRI4S:  Invalid digitization number!');
end
%
if (nargin<2)||isempty(kid)
  error(' *** ERROR in MK_TRI4S:  No knee identifier!');
end
%
if (nargin<1)||isempty(data)
  error(' *** ERROR in MK_TRI4S:  No input data!');
end
%
% Validate Knee Identifers
%
nc = size(kid,2);
if nc<3||nc>4
  error(' *** ERROR in MK_TRI4S:  Invalid knee identifier!');
end
if nc==3
  kid = [' ' kid];
end
%
vid = [' 4_R'; ' 5_L'; '56_R'; '57_L'; '64_R'; '66_R'];    % Valid knee identifers
%
idk = strmatch(upper(kid),vid);
if isempty(idk)||(idig==1&&idk==6)||(idig==2&&idk==1)||(idig==2&&idk==2)
  error(' *** ERROR in MK_TRI4S:  Invalid knee identifier!');
end
%
idk = idk+3*(idig-1);   % Unique number for each knee and digitization
%
% Check Knee Data and Get Knee Specific Indices
%
[nslice,nsl,isl] = sl_info(data);
%
if idk==1
  if nslice~=43||nsl(1)~=29||nsl(2)~=33||nsl(3)~=43||nsl(4)~=47|| ...
     isl(8)~=317
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  data = flipud(data);  % Reverse order of slices - only this knee
  [nslice,nsl,isl] = sl_info(data);
  id1 = 37;
elseif idk==2
  if nslice~=46||nsl(44)~=28||nsl(45)~=15||isl(nslice+1)~=3506
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 43;
elseif idk==3
  if nslice~=48||nsl(45)~=51||nsl(46)~=51||nsl(47)~=45|| ...
     isl(nslice+1)~=3841
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 43;
elseif idk==4
  if nslice~=41||nsl(39)~=53||nsl(40)~=19||isl(nslice+1)~=3212
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 37;
elseif idk==5
  if nslice~=46||nsl(43)~=62||nsl(44)~=57||nsl(45)~=47|| ...
     nsl(46)~=45||isl(nslice+1)~=4159
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 39;
elseif idk==6
  if nslice~=48
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 40;
elseif idk==7
  if nslice~=41||nsl(39)~=52||nsl(40)~=19||nsl(41)~=14|| ...
     isl(nslice+1)~=3170
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 36;
elseif idk==8
  if nslice~=46||nsl(44)~=59||nsl(45)~=48||isl(nslice)~=4113
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 40;
else
  if nslice~=40||nsl(39)~=51||nsl(40)~=43||isl(nslice+1)~=3424
    error(' *** ERROR in MK_TRI4S:  Invalid data!');
  end
  id1 = 37;
end
%
id2 = (id1:nslice)';
id1 = (1:id1)';
%
% Mesh Axial Data in Two Parts
%
dat1 = data(id1);
dat2 = data(id2);
%
% Mesh First Part
%
[tri1,xyz1] = mk_tri4p2(dat1);
%
nid = size(xyz1,1);     % Number of points in part1
nt1 = size(tri1,1);     % Number of triangles in part1
isl = isl(id1(end))+1:isl(id1(end)+1); % Nodes in bottom slice
%
% Remove Bottom Triangles and Bottom Node
%
it1 = nod2tri(nid,tri1);               % Triangles attached to bottom node
it2 = nod2tri(isl,tri1,2);             % Triangles attached only to bottom slice
it = unique([it1;it2]);
idx = true(nt1,1);
idx(it) = false;
tri1 = tri1(idx,:);     % Remove bottom triangles
%
nid = nid-1;            % Remove bottom node
xyz1 = xyz1(1:nid,:);   % Remove bottom node coordinates
%
% Mesh Second Part
%
tri2 = mk_tri4i(dat2,true,true,2);
%
% Correct Mesh
%
if idk==1
%
  tri2([137; 185],:) = [ 63    62   165
                        165   118    63 ];       % Corrected connectivity
%
  it = [283 284 289 240:244]';         % Triangles needing fixing
%
  tri2(it,:) = [167   164   163
                163   166   167
                164   168   169
                169   165   164
                165   169   170
                164   167   168
                175   176   125
                170   118   165 ];     % Corrected connectivity
%
  it = (317:544)';      % Triangles between slices 40 and 43
%
  tri2(it,:) = [227   186   226
                202   239   240
                219   179   178
                201   200   238
                214   175   174
                178   177   218
                182   222   223
                183   223   224
                237   199   198
                205   204   241
                189   230   190
                231   192   191
                207   243   244
                226   186   185
                167   166   251
                171   170   254
                172   171   255
                255   173   172
                168   167   252
                253   168   252
                253   169   168
                173   255   213
                219   178   218
                181   180   221
                222   181   221
                185   225   226
                216   177   176
                213   174   173
                224   184   183
                227   187   186
                228   188   187
                230   191   190
                196   195   235
                197   196   235
                239   202   201
                203   240   241
                206   242   243
                244   208   207
                212   248   249
                166   249   250
                241   242   205
                242   206   205
                208   244   245
                209   245   246
                210   246   247
                247   211   210
                200   199   237
                198   197   236
                232   193   192
                228   229   188
                234   195   194
                231   191   230
                194   233   234
                236   237   198
                238   239   201
                237   238   200
                240   203   202
                204   203   241
                211   247   248
                251   166   250
                245   209   208
                246   210   209
                248   212   211
                249   166   212
                243   207   206
                235   236   197
                195   234   235
                193   232   233
                232   192   231
                233   194   193
                220   179   219
                225   185   184
                229   189   188
                230   189   229
                228   187   227
                182   181   222
                180   220   221
                180   179   220
                184   224   225
                223   183   182
                171   254   255
                169   253   254
                174   213   214
                176   215   216
                175   214   215
                177   217   218
                177   216   217
                215   176   175
                167   251   252
                254   170   169
                249   248   277
                253   252   280
                240   269   270
                238   237   268
                262   228   261
                224   223   258
                222   256   257
                214   282   283
                254   281   255
                250   278   279
                281   254   253
                280   252   251
                275   246   274
                276   248   247
                279   251   250
                278   250   249
                277   248   276
                276   247   275
                269   239   238
                264   232   263
                267   237   236
                245   273   274
                272   242   271
                268   237   267
                226   261   227
                258   259   224
                256   222   221
                287   218   286
                257   223   222
                288   219   287
                221   220   288
                286   217   285
                216   284   285
                215   283   284
                282   214   213
                266   235   265
                264   234   233
                262   231   230
                260   225   259
                229   262   230
                263   231   262
                261   226   260
                219   288   220
                218   287   219
                217   286   218
                285   217   216
                284   216   215
                223   257   258
                288   256   221
                283   215   214
                213   281   282
                231   263   232
                232   264   233
                228   262   229
                261   228   227
                224   259   225
                225   260   226
                272   273   243
                269   240   239
                270   271   241
                265   235   234
                270   241   240
                271   242   241
                238   268   269
                235   266   236
                266   267   236
                234   264   265
                253   280   281
                278   249   277
                274   246   245
                275   247   246
                272   243   242
                213   255   281
                251   279   280
                273   244   243
                273   245   244
                289   257   290
                257   289   317
                266   265   310
                299   298   279
                286   285   294
                275   274   303
                294   293   286
                261   314   262
                314   261   260
                259   315   260
                260   315   314
                316   259   258
                315   259   316
                267   266   309
                312   263   313
                311   264   312
                310   265   311
                270   307   306
                269   308   307
                270   306   271
                272   306   305
                273   305   304
                273   304   274
                304   303   274
                276   275   302
                278   277   300
                282   281   280
                279   278   299
                284   296   295
                291   256   288
                287   286   293
                256   291   290
                293   292   287
                301   300   277
                277   276   301
                297   283   282
                296   284   283
                283   297   296
                298   297   280
                280   279   298
                280   297   282
                292   291   288
                288   287   292
                300   299   278
                295   294   285
                285   284   295
                273   272   305
                303   302   275
                302   301   276
                309   266   310
                308   269   268
                308   267   309
                272   271   306
                265   264   311
                268   267   308
                269   307   270
                290   257   256
                262   314   313
                264   263   312
                317   258   257
                258   317   316
                263   262   313 ];     % Corrected connectivity
%
elseif idk==2
  it = (61:103)';       % Triangles needing fixing
%
  tri2(it,:) = [ 73    58    72
                 58    73    59
                 72    58    57
                 57    71    72
                 71    57    56
                 56    70    71
                 70    56    55
                 69    70    55
                 69    55    54
                 68    69    54
                 68    54    53
                 68    53    52
                 67    68    52
                 67    52    51
                 66    67    51
                 66    51    50
                 65    66    50
                 65    50    49
                 65    49    48
                 48    64    65
                 64    48    47
                 47    63    64
                 63    47    46
                 63    46    45
                 45    62    63
                 62    45    44
                 44    36    62
                 36    44    43
                 73    74    59
                 74    60    59
                 74    75    60
                 75    33    60
                 33    61    34
                 61    33    75
                 61    35    34
                 61    62    35
                 62    36    35
                 36    43    37
                 37    42    38
                 38    41    39
                 39    41    40
                 38    42    41
                 37    43    42 ];     % Corrected connectivity
%
elseif idk==3
  it = [368 370 389 417 454 457 473]'; % Triangles needing fixing
%
  tri2(it,:) = [188   235   236
                234   235   187
                188   187   235
                277   236   235
                234   277   235
                276   234   233
                277   234   276 ];     % Corrected connectivity
%
elseif idk==4
  it = (218:289)';      % Triangles between slices 39 and 40
%
  tri2(it,:) = [120   153   152
                154   118   117
                119   154   153
                155   117   116
                115   157   156
                158   113   159
                113   158   114
                125   150   168
                112   159   113
                151   122   121
                121   152   151
                122   151   123
                152   121   120
                123   151   124
                154   119   118
                153   120   119
                116   156   155
                117   155   154
                163   160   159
                156   116   115
                163   162   161
                160   163   161
                159   112   111
                114   158   157
                157   115   114
                124   151   150
                150   125   124
                126   168   127
                176   136   135
                137   136   177
                129   128   170
                175   133   174
                135   134   175
                181   182   143
                143   142   181
                175   176   135
                144   143   182
                145   144   182
                166   148   147
                130   171   172
                133   175   134
                139   178   179
                182   164   145
                164   146   145
                136   176   177
                177   178   137
                180   181   142
                142   141   180
                140   179   180
                141   140   180
                146   164   165
                147   146   165
                165   166   147
                174   132   173
                178   138   137
                178   139   138
                179   140   139
                172   131   130
                173   131   172
                129   171   130
                132   174   133
                131   173   132
                170   171   129
                128   169   170
                169   128   127
                127   168   169
                168   126   125
                168   150   149
                149   167   168
                167   149   148
                148   166   167
                111   163   159 ];
%
elseif idk==5
%
  it = [200 236 244 248 263 264];      % Triangles needing fixing
  tri2(it,:) = [134   200   201
                134   133   200
                133   132   199
                133   199   200
                131   198   199
                132   131   199];      % Corrected connectivity
%
  it = [552 569 571 577 610 628 634]'; % Triangles needing fixing 
  tri2(it,:) = [352   347   351
                351   347   350
                350   347   346
                347   352   286
                345   348   349
                349   346   345
                346   349   350 ];     % Corrected connectivity
%
  it = [681 685 689 690 692 693 715 716 717 718 719 720 721 722 723 ...
        724 725 737 738 744 767 768 769 772 773 774 780 789 794 795 ... 
        796 797 798 799 800 801 802 807 823 832 833 834 835 836 837 ...
        841]';          % Triangles needing fixing 
%
  tri2(it,:) = [391   405   390
                391   405   392
                392   406   393
                393   407   394
                394   408   395
                395   408   396
                396   410   397
                397   411   398
                398   411   399
                399   412   400
                400   413   401
                406   392   405
                407   393   406
                408   394   407
                396   408   409
                410   396   409
                411   397   410
                399   411   412
                400   412   413
                401   413   414
                401   414   402
                402   415   403
                403   416   404
                415   402   414
                404   416   348
                416   403   415
                405   489   406
                406   490   407
                407   491   408
                408   492   409
                409   493   410
                410   494   411
                411   495   412
                413   495   496
                405   488   489
                406   489   490
                407   490   491
                408   491   492
                409   492   493
                410   493   494
                411   494   495
                495   413   412
                413   496   414
                414   496   415
                415   452   416
                415   496   452 ];     % Corrected connectivity
%
  it = [614 618 678 730 740 746];      % Triangles needing fixing 
  tri2(it,:) = [386   325   385
                325   384   385
                386   447   387
                386   385   446
                447   386   446
                384   446   385 ];     % Corrected connectivity
%
  it = [775 777 814 825 829 838 842];
  tri2(it,:) = [446   445   480
                445   479   480
                445   444   479
                444   478   479
                444   443   478
                443   477   478
                443   442   477];
%
elseif idk==7
  it = (326:396)';      % Triangles between slices 39 and 40
%
  tri2(it,:) = [210   177   211
                212   175   174
                177   176   211
                174   213   212
                171   215   172
                215   171   170
                170   169   216
                207   226   183
                217   169   168
                208   207   182
                209   178   210
                208   179   209
                179   178   209
                208   180   179
                176   175   211
                178   177   210
                214   173   172
                212   211   175
                217   168   218
                213   174   173
                168   167   218
                216   169   217
                216   215   170
                215   214   172
                214   213   173
                181   180   208
                208   182   181
                182   207   183
                191   232   233
                192   234   193
                184   227   185
                231   189   230
                232   191   190
                237   219   199
                199   198   237
                190   231   232
                200   199   219
                220   201   200
                203   222   223
                186   227   228
                189   231   190
                235   195   194
                200   219   220
                201   220   221
                233   192   191
                234   192   233
                236   237   196
                198   197   237
                236   196   195
                196   237   197
                221   202   201
                222   203   202
                202   221   222
                230   188   229
                234   194   193
                235   194   234
                236   195   235
                228   187   186
                229   187   228
                185   227   186
                188   230   189
                187   229   188
                226   227   184
                226   184   183
                226   207   206
                225   226   206
                225   206   205
                224   225   205
                224   205   204
                223   224   204
                223   204   203 ];     % Corrected connectivity
%
elseif idk==8
  it = [535 550 607]';  % Triangles needing fixing
  tri2(it,:) = [289   343   290
                343   291   290
                291   343   344 ];     % Corrected connectivity
%
elseif idk==9
  it = [230 236 268 281 283 287]';     % Triangles needing fixing
  tri2(it,:) = [154   198   155
                156   155   198
                197   154   153
                153   196   197
                198   154   197
                198   199   156 ];     % Corrected connectivity
end
%
% Get Coordinates
%
xyz2 = cell2mat(dat2);
%
if iplt
  figure;
  orient tall;
  tri_shr(tri2,xyz2);
  view(3);
  axis equal;
end
%
% Combine Meshes
%
nn = nsl(id2(1));
idn = find(tri2>nn);
tri2(idn) = tri2(idn)+1;               % Offset by one for addition of top node
nn = nn+1;
tria = [tri1; tri2+(nid-nn)];
xyza = [xyz1; xyz2(nn:end,:)];
%
if nargout>3
  nt = size(tria,1);
end
%
if iplt
  figure;
  orient tall;
  tri_shr(tria,xyza);
  view(3);
  axis equal;
end
%
return