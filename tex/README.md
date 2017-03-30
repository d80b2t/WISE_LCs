# README for d80b2t/WISE_LCs/tex

Okay, I (NPR) haven't quite got my full act together, and the .tex and subsequent write-up
of the WISE QSOs and Light Curve project(s) isn't fully on this github repo (yet).

However, all the relevant documentation and files can currently be found on Overleaf:
https://www.overleaf.com/7362842xzrqrnwmxsyz#/25578795/



## emulateapj.cls vs. revtex4-1 vs. lscape

tl; dr: Disable \usepackage{lscape}


On May 3, 2011, at 7:09 AM, Jake Weiskoff [www-admin] <www-admin@arxiv.org> wrote:

Dear Nic,

This is because of a known incompatibility with lscape/rotating and revtex4-1 which is used by emulateapj.cls. There are numerous discussions of this in various tex forums, but the simple fix is to force your document to use revtex4. I've done this by updating your class file:

0240675/src> diff -u emulateapj.cls emulateapj.old
--- emulateapj.cls	2011-05-03 10:05:55.759528050 -0400
+++ emulateapj.old	2011-05-02 17:39:09.000000000 -0400
@@ -70,8 +70,7 @@
\newif\if@revtex@four@one@
\IfFileExists{revtex4-1.cls}{
\@revtex@four@one@true
-%\def\@revtex@cls{revtex4-1}
-\def\@revtex@cls{revtex4}
+\def\@revtex@cls{revtex4-1}
}{
\@revtex@four@one@false
\def\@revtex@cls{revtex4}


Please clear your browser's cache and check your submission.

--
arXiv admin



	Wed Mar 29 11:39:51 PDT 2017

