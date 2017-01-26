;;; Spheral-style.el
;;;
;;; Based on kull-style.el from Kat Price

;;; A standard style for Spheral.

(defun indent-buffer ()
  "Indent the entire buffer"
  (interactive "")
  (indent-region (point-min) (point-max) nil))

(defun spheral-c-mode-common-hook ()
  ;; Use gnu style but with no extra substatement braces indentation 
  (c-set-style "gnu")
  (setq c-basic-offset 2)
  (c-set-offset 'substatement-open 0)
  ;; Set width of a tab character to be 6 spaces
  ;;  (setq tab-width 6)
  ;; Inserting a tab really inserts spaces instead of a tab
  (setq indent-tabs-mode nil)
  ;; Setup return key to insert a newline and indent it
  (define-key c-mode-map "\C-m" 'newline-and-indent)
  )
(add-hook 'c-mode-common-hook 'spheral-c-mode-common-hook)
(remove-hook 'c-mode-common-hook 'kull-c-mode-common-hook)

;;; Provide the standard function intro separator
(defun spheral-insert-function-doc ()
  (interactive)
  (beginning-of-line)
  (insert "//------------------------------------------------------------------------------\n")
  (insert "// \n")
  (insert "//------------------------------------------------------------------------------\n\n")
  (previous-line 3)
  (end-of-line)
)

;;; Make Python code use 4 space indentation.  Also, force Emacs to
;;; use spaces instead of tabs, since mixed mode causes problems with
;;; some editors.

(defun spheral-python-mode-hook ()
  (setq py-indent-offset 4)
  ;; Inserting a tab really inserts spaces instead of a tab
  (setq indent-tabs-mode nil)
  ;; Setup return key to insert a newline and indent it
  (define-key c-mode-map "\C-m" 'newline-and-indent)
  )
(add-hook 'python-mode-hook 'spheral-python-mode-hook)
(remove-hook 'python-mode-common-hook 'kull-python-mode-common-hook)

;;; Define some standard key shortcuts.
(global-set-key [f1] 'compile)
(global-set-key [f7] 'spheral-insert-function-doc)
