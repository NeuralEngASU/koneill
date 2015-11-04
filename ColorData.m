%%%%
% Adapted from research:
% ColorBrewer.org: An Online Tool for Selecting Colour
% Schemes for Maps
% Brewer, C
%
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.361.6082&rep=rep1&type=pdf
%
%
% 5. Citation
% Ahhh, that lavish citation I'm hoping for...
% Wording depends on context and I'm not too picky. If you add a note in the corner of a map, how about one of these:
%     - Colors from www.ColorBrewer.org by Cynthia A. Brewer, Geography, Pennsylvania State University.
%     - Map colors based on www.ColorBrewer.org, by Cynthia A. Brewer, Penn State.
%     - Color symbols: ColorBrewer.org
% A reference in a journal article or report might look like:
%     - Brewer, Cynthia A., 200x. http://www.ColorBrewer.org, accessed date.
%%%

%% Sequential Colors
YlGn = [255	255	229
    247	252	185
    217	240	163
    173	221	142
    120	198	121
    65	171	93
    35	132	67
    0	104	55
    0	69	41
    ];

YlGnBu = [255	255	217
    237	248	177
    199	233	180
    127	205	187
    65	182	196
    29	145	192
    34	94	168
    37	52	148
    8	29	88
    ];

BuGn = [247	252	240
224	243	219
204	235	197
168	221	181
123	204	196
78	179	211
43	140	190
8	104	172
8	64	129];

GnBu = [247	252	253
229	245	249
204	236	230
153	216	201
102	194	164
65	174	118
35	139	69
0	109	44
0	68	27];

PuBuGn = [255	247	251
236	226	240
208	209	230
166	189	219
103	169	207
54	144	192
2	129	138
1	108	89
1	70	54];

PuBu = [255	247	251
236	231	242
208	209	230
166	189	219
116	169	207
54	144	192
5	112	176
4	90	141
2	56	88];

BuPu = [247	252	253
224	236	244
191	211	230
158	188	218
140	150	198
140	107	177
136	65	157
129	15	124
77	0	75];

RdPu = [255	247	243
253	224	221
252	197	192
250	159	181
247	104	161
221	52	151
174	1	126
122	1	119
73	0	106];

PuRd = [247	244	249
231	225	239
212	185	218
201	148	199
223	101	176
231	41	138
206	18	86
152	0	67
103	0	31];

OrRd = [255	247	236
254	232	200
253	212	158
253	187	132
252	141	89
239	101	72
215	48	31
179	0	0
127	0	0];

YlOrRd = [255	255	204
255	237	160
254	217	118
254	178	76
253	141	60
252	78	42
227	26	28
189	0	38
128	0	38];

YlOrBr = [255	255	229
255	247	188
254	227	145
254	196	79
254	153	41
236	112	20
204	76	2
153	52	4
102	37	6];

Pu = [252	251	253
239	237	245
218	218	235
188	189	220
158	154	200
128	125	186
106	81	163
84	39	143
63	0	125];

Bu = [247	251	255
222	235	247
198	219	239
158	202	225
107	174	214
66	146	198
33	113	181
8	81	156
8	48	107];

Gn = [247	252	245
229	245	224
199	233	192
161	217	155
116	196	118
65	171	93
35	139	69
0	109	44
0	68	27];

Or = [255	245	235
254	230	206
253	208	162
253	174	107
253	141	60
241	105	19
217	72	1
166	54	3
127	39	4];

Rd = [255	245	240
254	224	210
252	187	161
252	146	114
251	106	74
239	59	44
203	24	29
165	15	21
103	0	13];



%%
YlGn = YlGn./255;

%%

hold on
plot(YlGn(:,1), 'r')
plot(YlGn(:,2), 'g')
plot(YlGn(:,3), 'b')
hold off


%%


len = 11;
n = 5;

p1 = polyfit([1:len]', YlGn(:,1), n);
p2 = polyfit([1:len]', YlGn(:,2), n);
p3 = polyfit([1:len]', YlGn(:,3), n);

val1 = zeros(100,3);

temp1 = p1(1) + p1(2) + p1(3) + p1(4);
temp2 = p2(1) + p2(2) + p2(3) + p2(4);
temp3 = p3(1) + p3(2) + p3(3) + p3(4);

for i=1:100
    for j = 1:n+1
    
    k = n:-1:0;
    x = linspace(1,len,100);
    val1(i,1) = val1(i,1) + p1(j)*x(i)^(k(j));
    val1(i,2) = val1(i,2) + p2(j)*x(i)^(k(j));
    val1(i,3) = val1(i,3) + p3(j)*x(i)^(k(j));
    

    end % END FOR
    
    if val1(i, 1) > 1
        val1(i,1) = 1;
    end
    if val1(i, 2) > 1
        val1(i,2) = 1;
    end
    if val1(i, 3) > 1
        val1(i,3) = 1;
    end
    
    if val1(i, 1) < 0
        val1(i,1) = 0;
    end
    if val1(i, 2) < 0
        val1(i,2) = 0;
    end
    if val1(i, 3) < 0
        val1(i,3) = 0;
    end
    
end %END FOR


hold on
plot(x,val1(:,1), 'r')
plot(x,val1(:,2), 'g')
plot(x,val1(:,3), 'b')

plot(YlGn(:,1), 'r:')
plot(YlGn(:,2), 'g:')
plot(YlGn(:,3), 'b:')
hold off


%%

img = zeros(100,50,3);


for i = 1:100
    
    
    img(i, :, 1) = repmat(val1(i, 1), 1,50);
    img(i, :, 2) = repmat(val1(i, 2), 1,50);
    img(i, :, 3) = repmat(val1(i, 3), 1,50);
    
    
end % END FOR

imshow(img)

figure (2)

ttemp = 1:0.01:10;
wtemp = rand(1, length(ttemp));
xtemp = sin(ttemp);
ytemp = cos(ttemp);
ztemp = exp(ttemp./5);

hold on
ldataw = plot(ttemp, wtemp);
ldatax = plot(ttemp, xtemp);
ldatay = plot(ttemp, ytemp);
ldataz = plot(ttemp, ztemp);
hold off
set(ldataw, 'LineWidth', 2.75);
set(ldatax, 'LineWidth', 2.75);
set(ldatay, 'LineWidth', 2.75);
set(ldataz, 'LineWidth', 2.75);


set(ldataw, 'Color', val1(25,:))
set(ldatax, 'Color', val1(15,:))
set(ldatay, 'Color', val1(75,:))
set(ldataz, 'Color', val1(80,:))

%% Diverging Colors

PuOr = [127	59	8
179	88	6
224	130	20
253	184	99
254	224	182
247	247	247
216	218	235
178	171	210
128	115	172
84	39	136
45	0	75];

BrGn = [64	0	75
118	42	131
153	112	171
194	165	207
231	212	232
247	247	247
217	240	211
166	219	160
90	174	97
27	120	55
0	68	27];

BuOr = [165	0	38
215	48	39
244	109	67
253	174	97
254	224	144
255	255	191
224	243	248
171	217	233
116	173	209
69	117	180
49	54	149];

PRGn = [64	0	75
118	42	131
153	112	171
194	165	207
231	212	232
247	247	247
217	240	211
166	219	160
90	174	97
27	120	55
0	68	27];

%% Qualaitative

numPoints = 3; % 3-12

% Sequence 1
% qual1 = [ 84,  48,   5;...
%     140,  81,  10;...
%     191, 129,  45;...
%     223, 194, 125;...
%     246, 232, 195;...
%     245, 245, 245;...
%     199, 234, 229;...
%     128, 205, 193;...
%     53, 151, 143;...
%     1, 102,  94;...
%     0,  60,  48
    
set1 = [141	211	199
255	255	179
190	186	218
251	128	114
128	177	211
253	180	98
179	222	105
252	205	229
217	217	217
188	128	189
204	235	197
255	237	111];

% qual1 =
% #a6cee3
% #1f78b4
% #b2df8a
% #33a02c
% #fb9a99
% #e31a1c
% #fdbf6f
% #ff7f00
% #cab2d6
% #6a3d9a
% #ffff99


% qual2 =[141, 211, 199;...
%     255, 255, 179;...
%     190, 186, 218;...
%     251, 128, 114;...
%     128, 177, 211;...
%     253, 180,  98;...
%     179, 222, 105;...
%     252, 205, 229;...
%     217, 217, 217;...
%     188, 128, 189;...
%     204, 235, 197];

% qual2 =
% #8dd3c7
% #ffffb3
% #bebada
% #fb8072
% #80b1d3
% #fdb462
% #b3de69
% #fccde5
% #d9d9d9
% #bc80bd
% #ccebc5