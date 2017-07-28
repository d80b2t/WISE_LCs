%% Class nature
%% Written by Peter Czoschke, czoschke@mrl.uiuc.edu
%% v1.0 24 Feb. 2004
%%
%% A document class for Nature letters and articles that (should)
%% conform to the Nature preprint style as defined by their Word template.
%%
%% I created this class for personal purposes and without any connection
%% to the Nature Publishing Group.  I in no way claim that documents generated
%% with this file fully comply with their current style requirements.
%% I disclaim any responsibility for the use of this file heretofore.
%%
%% ------------------------------------
%%
%% See the accompanying file nature-template.tex for an example.
%%
%% Use with naturemag.bst if using BibTeX.
%%
%% Nature wants things in the order:
%% title, text, methods, references, Supplementary Information line (if any),
%% acknowledgements, interest declaration, corresponding author line,
%% tables, figure legends
%%
%% This class provides the following non-standard features:
%% 
%%  1   "affiliations" environment for listing the authors' institutions
%%
%%  2   "methods" environment for the Methods section
%%
%%  3   "addendum" environment for the Acknowledgements, etc. at the end
%%
%% In addition, the abstract, figure and table environments are redefined. See below.
%%

\ProvidesClass{nature}[24/2/2004 v1.0]
\typeout{A class for writing preprints for the journal Nature}
\typeout{Written by Peter Czoschke}
\typeout{ }

\LoadClass[10pt]{article}
\RequirePackage{cite}
\RequirePackage{times}
\RequirePackage{fullpage}
\RequirePackage{ifthen}

%% make labels in bibliobraphy be #.
\renewcommand\@biblabel[1]{#1.}

%% make citations be superscripts, taken from citesupernumber.sty
\def\@cite#1#2{$^{\mbox{\scriptsize #1\if@tempswa , #2\fi}}$}

%% Some style parameters
\setlength{\parindent}{0.39in}
\setlength{\parskip}{18pt}
%\newcommand{\spacing}[1]{\renewcommand{\baselinestretch}{#1}\large\normalsize}
\newcommand{\spacing}[1]{\renewcommand{\baselinestretch}{1.2}\large\normalsize}
\spacing{2}

%% Redefine \maketitle for Nature style
\def\@maketitle{%
  \newpage\spacing{1}\setlength{\parskip}{12pt}%
    {\Large\bfseries\noindent\sloppy \textsf{\@title} \par}%
    {\noindent\sloppy \@author}%
}

%% Define the affiliations environment, list each institution as an \item
%% Put after \maketitle
\newenvironment{affiliations}{%
    \setcounter{enumi}{1}%
    \setlength{\parindent}{0in}%
    \slshape\sloppy%
    \begin{list}{\upshape$^{\arabic{enumi}}$}{%
        \usecounter{enumi}%
        \setlength{\leftmargin}{0in}%
        \setlength{\topsep}{0in}%
        \setlength{\labelsep}{0in}%
        \setlength{\labelwidth}{0in}%
        \setlength{\listparindent}{0in}%
        \setlength{\itemsep}{0ex}%
        \setlength{\parsep}{0in}%
        }
    }{\end{list}\par\vspace{12pt}}

%% Redefine the abstract environment to be the first bold paragraph
\renewenvironment{abstract}{%
    \setlength{\parindent}{0in}%
    \setlength{\parskip}{0in}%
    \bfseries%
    }{\par\vspace{-6pt}}

%% Redefine the \section command.
\renewcommand{\section}{\@startsection {section}{1}{0pt}%
    {-6pt}{1pt}%
    {\bfseries\centering}%
    }
\renewcommand{\subsection}{\@startsection {subsection}{2}{0pt}%
    {-0pt}{-0.5em}%
    {\bfseries}*%
    }

%% Define the methods environment.  Use \subsection to separate.
\newenvironment{methods}{%
    \section*{Methods}%
    \setlength{\parskip}{12pt}%
    }{}

%% No heading for References section, but eat up the extra space from \section command
\renewcommand\refname{\vspace{-48pt}\setlength{\parskip}{12pt}}

%% Define the addendum environment for Supplementary Info, Acknowledgements, etc.
\newenvironment{addendum}{%
    \setlength{\parindent}{0in}%
    \small%
    \begin{list}{Acknowledgements}{%
        \setlength{\leftmargin}{0in}%
        \setlength{\listparindent}{0in}%
        \setlength{\labelsep}{0em}%
        \setlength{\labelwidth}{0in}%
        \setlength{\itemsep}{12pt}%
        \let\makelabel\addendumlabel}
    }
    {\end{list}\normalsize}

\newcommand*{\addendumlabel}[1]{\textbf{#1}\hspace{1em}}

%% Redefine the figure and table environments.
%% The new figure environment defers printing of all figure captions
%% until the end of the document and prints them in sans serif font.
%% The \includegraphics command is ignored since Nature wants figures
%% as separate files.
%% The table environment was just changed to print on a separate page in
%% sans serif font.  Tables should all be moved to the end of the document
%% manually.  

\usepackage[labelsep=space,labelfont=bf,textfont=md]{caption}


% \let\svfigure\figure
%\let\endsvfigure\endfigure

%\renewenvironment{figure}{\let\caption\NAT@figcaption\svfigure}{\endsvfigure}
%\renewenvironment{figure*}{\let\caption\NAT@figcaption\svfigure}{\endsvfigure}
 %
%\renewenvironment{figure}{\let\caption\NAT@figcaption}{}
%\renewenvironment{figure*}{\let\caption\fullwidthcaption}{}
    
\newcommand{\NAT@figcaption}[2][]{%
\onecolumn  
    \refstepcounter{figure}
    \ifthenelse{\value{figure}=1}{
        \newpage\noindent%
        %\rule{\textwidth}{1pt}
    }{
        \par\vfill
    }
    \sffamily\noindent\textbf{Figure \arabic{figure}}\hspace{0.1em}#2
  
    }

\newcommand{\fullwidthcaption}[2][]{%
    \refstepcounter{figure}
    \ifthenelse{\value{figure}=1}{
        \newpage\noindent%
        %\rule{\textwidth}{1pt}
    }{
        \par\vfill
    }
    \sffamily\noindent\textbf{Figure \arabic{figure}}\hspace{0.1em}#2
    }

\renewenvironment{table}{\@float{table}[p]\sffamily}{\end@float}

%% ignore all \includegraphics commands in the document
%\AtBeginDocument{\let\includegraphics\NAT@ignore}
%\newcommand{\NAT@ignore}[2][]{}

\let\jnl@style=\rm
\def\ref@jnl#1{{\jnl@style#1}}

\def\aj{\ref@jnl{AJ}}                   % Astronomical Journal
\def\actaa{\ref@jnl{Acta Astron.}}      % Acta Astronomica
\def\araa{\ref@jnl{ARA\&A}}             % Annual Review of Astron and Astrophys
\def\apj{\ref@jnl{ApJ}}                 % Astrophysical Journal
\def\apjl{\ref@jnl{ApJ}}                % Astrophysical Journal, Letters
\def\apjs{\ref@jnl{ApJS}}               % Astrophysical Journal, Supplement
\def\ao{\ref@jnl{Appl.~Opt.}}           % Applied Optics
\def\apss{\ref@jnl{Ap\&SS}}             % Astrophysics and Space Science
\def\aap{\ref@jnl{A\&A}}                % Astronomy and Astrophysics
\def\aapr{\ref@jnl{A\&A~Rev.}}          % Astronomy and Astrophysics Reviews
\def\aaps{\ref@jnl{A\&AS}}              % Astronomy and Astrophysics, Supplement
\def\azh{\ref@jnl{AZh}}                 % Astronomicheskii Zhurnal
\def\baas{\ref@jnl{BAAS}}               % Bulletin of the AAS
\def\bac{\ref@jnl{Bull. astr. Inst. Czechosl.}}
                % Bulletin of the Astronomical Institutes of Czechoslovakia 
\def\caa{\ref@jnl{Chinese Astron. Astrophys.}}
                % Chinese Astronomy and Astrophysics
\def\cjaa{\ref@jnl{Chinese J. Astron. Astrophys.}}
                % Chinese Journal of Astronomy and Astrophysics
\def\icarus{\ref@jnl{Icarus}}           % Icarus
\def\jcap{\ref@jnl{J. Cosmology Astropart. Phys.}}
                % Journal of Cosmology and Astroparticle Physics
\def\jrasc{\ref@jnl{JRASC}}             % Journal of the RAS of Canada
\def\memras{\ref@jnl{MmRAS}}            % Memoirs of the RAS
\def\mnras{\ref@jnl{MNRAS}}             % Monthly Notices of the RAS
\def\na{\ref@jnl{New A}}                % New Astronomy
\def\nar{\ref@jnl{New A Rev.}}          % New Astronomy Review
\def\pra{\ref@jnl{Phys.~Rev.~A}}        % Physical Review A: General Physics
\def\prb{\ref@jnl{Phys.~Rev.~B}}        % Physical Review B: Solid State
\def\prc{\ref@jnl{Phys.~Rev.~C}}        % Physical Review C
\def\prd{\ref@jnl{Phys.~Rev.~D}}        % Physical Review D
\def\pre{\ref@jnl{Phys.~Rev.~E}}        % Physical Review E
\def\prl{\ref@jnl{Phys.~Rev.~Lett.}}    % Physical Review Letters
\def\pasa{\ref@jnl{PASA}}               % Publications of the Astron. Soc. of Australia
\def\pasp{\ref@jnl{PASP}}               % Publications of the ASP
\def\pasj{\ref@jnl{PASJ}}               % Publications of the ASJ
\def\rmxaa{\ref@jnl{Rev. Mexicana Astron. Astrofis.}}%
                % Revista Mexicana de Astronomia y Astrofisica
\def\qjras{\ref@jnl{QJRAS}}             % Quarterly Journal of the RAS
\def\skytel{\ref@jnl{S\&T}}             % Sky and Telescope
\def\solphys{\ref@jnl{Sol.~Phys.}}      % Solar Physics
\def\sovast{\ref@jnl{Soviet~Ast.}}      % Soviet Astronomy
\def\ssr{\ref@jnl{Space~Sci.~Rev.}}     % Space Science Reviews
\def\zap{\ref@jnl{ZAp}}                 % Zeitschrift fuer Astrophysik
\def\nat{\ref@jnl{Nature}}              % Nature
\def\iaucirc{\ref@jnl{IAU~Circ.}}       % IAU Cirulars
\def\aplett{\ref@jnl{Astrophys.~Lett.}} % Astrophysics Letters
\def\apspr{\ref@jnl{Astrophys.~Space~Phys.~Res.}}
                % Astrophysics Space Physics Research
\def\bain{\ref@jnl{Bull.~Astron.~Inst.~Netherlands}} 
                % Bulletin Astronomical Institute of the Netherlands
\def\fcp{\ref@jnl{Fund.~Cosmic~Phys.}}  % Fundamental Cosmic Physics
\def\gca{\ref@jnl{Geochim.~Cosmochim.~Acta}}   % Geochimica Cosmochimica Acta
\def\grl{\ref@jnl{Geophys.~Res.~Lett.}} % Geophysics Research Letters
\def\jcp{\ref@jnl{J.~Chem.~Phys.}}      % Journal of Chemical Physics
\def\jgr{\ref@jnl{J.~Geophys.~Res.}}    % Journal of Geophysics Research
\def\jqsrt{\ref@jnl{J.~Quant.~Spec.~Radiat.~Transf.}}
                % Journal of Quantitiative Spectroscopy and Radiative Transfer
\def\memsai{\ref@jnl{Mem.~Soc.~Astron.~Italiana}}
                % Mem. Societa Astronomica Italiana
\def\nphysa{\ref@jnl{Nucl.~Phys.~A}}   % Nuclear Physics A
\def\physrep{\ref@jnl{Phys.~Rep.}}   % Physics Reports
\def\physscr{\ref@jnl{Phys.~Scr}}   % Physica Scripta
\def\planss{\ref@jnl{Planet.~Space~Sci.}}   % Planetary Space Science
\def\procspie{\ref@jnl{Proc.~SPIE}}   % Proceedings of the SPIE

\let\astap=\aap
\let\apjlett=\apjl
\let\apjsupp=\apjs
\let\applopt=\ao

\def\micron{\mathrel{\mu\rm m}} 
\def\arcsec{"} 
\def\lt{\mathrel{<}} 
\def\gt{\mathrel{>}} 
\def\tilde{{$\sim$}} 
\def\gsim{{$\gtrsim$}} 
\def\sime{{$\sim$}} 
\makeatletter
\newcommand{\ion}[2]{#1$\;$\textsc{\rmfamily\@roman{#2}}\relax}
\makeatother