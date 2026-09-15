#!/bin/bash
# Load full shell profile so npm/nvm are available
export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && source "$NVM_DIR/nvm.sh"
[ -f "$HOME/.bashrc" ] && source "$HOME/.bashrc" 2>/dev/null

echo "=============================="
echo "  AmpliconApp Frontend :5173"
echo "=============================="
echo ""
cd ~/r16s-app/frontend
npm run dev
echo ""
echo "Frontend stopped. Press Enter to close..."
read
