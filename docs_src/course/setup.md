# Setup Guide

## Before You Arrive

Complete these setup steps **at least one week** before the course begins.

## 1. System Requirements

### Minimum Specifications
- **RAM**: 8GB (16GB recommended)
- **Storage**: 50GB free space
- **OS**: Linux, macOS, or Windows 10/11
- **Internet**: Stable broadband connection

### Operating System Setup

=== "Linux (Ubuntu/Debian)"
    ```bash
    # Update package list
    sudo apt update && sudo apt upgrade -y
    
    # Install essential tools
    sudo apt install -y git curl wget build-essential
    ```

=== "macOS"
    ```bash
    # Install Homebrew (if not already installed)
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    
    # Install essential tools
    brew install git curl wget
    ```

=== "Windows"
    1. **Install Git Bash**: Download from [git-scm.com](https://git-scm.com/)
    2. **Enable WSL2** (optional but recommended):
       ```powershell
       wsl --install
       ```
    3. **Install Windows Terminal** from Microsoft Store

## 2. Git Configuration

Configure Git with your information:

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
git config --global init.defaultBranch main
```

Verify configuration:
```bash
git config --list
```

## 3. SSH Setup

Generate SSH key for secure connections:

```bash
# Generate SSH key pair
ssh-keygen -t ed25519 -C "your.email@example.com"

# Start SSH agent
eval "$(ssh-agent -s)"

# Add key to agent
ssh-add ~/.ssh/id_ed25519

# Display public key (copy this for later use)
cat ~/.ssh/id_ed25519.pub
```

## 4. Text Editor

Choose and install a text editor:

=== "VS Code (Recommended)"
    - Download from [code.visualstudio.com](https://code.visualstudio.com/)
    - Install useful extensions:
      - Python
      - Markdown All in One
      - GitLens

=== "Sublime Text"
    - Download from [sublimetext.com](https://www.sublimetext.com/)
    - Consider Sublime Merge for Git integration

=== "Command Line Editors"
    ```bash
    # For vim users
    sudo apt install vim  # Linux
    brew install vim      # macOS
    
    # For nano users (usually pre-installed)
    which nano
    ```

## 5. Test Your Setup

### Basic Command Line Test
```bash
# Check versions
git --version
curl --version
wget --version

# Test directory operations
mkdir test_setup
cd test_setup
echo "Hello World" > test.txt
cat test.txt
cd ..
rm -rf test_setup
```

### Git Test
```bash
# Clone a test repository
git clone https://github.com/CIDRI-Africa/microbial-genomics-training.git
cd microbial-genomics-training
ls -la
cd ..
rm -rf microbial-genomics-training
```

## 6. Course-Specific Setup

### HPC Access (If Provided)
You will receive:
- **SSH connection details**
- **Username and login instructions**
- **VPN setup** (if required)

### Docker/Singularity (Optional)
For local analysis (advanced users):

=== "Docker"
    ```bash
    # Linux
    sudo apt install docker.io
    sudo usermod -aG docker $USER
    
    # macOS
    # Download Docker Desktop from docker.com
    
    # Test installation
    docker --version
    docker run hello-world
    ```

=== "Singularity"
    ```bash
    # Linux (Ubuntu/Debian)
    sudo apt install singularity-container
    
    # Test installation
    singularity --version
    ```

## 7. Download Course Materials

One week before the course:

```bash
# Create course directory
mkdir -p ~/microbial-genomics-course
cd ~/microbial-genomics-course

# Clone course repository
git clone https://github.com/CIDRI-Africa/microbial-genomics-training.git

# Check repository contents
cd microbial-genomics-training
ls -la
```

## 8. Pre-course Data Download

Large datasets will be provided via:
- **Shared storage** on HPC systems
- **Cloud storage links** (Google Drive/Dropbox)
- **Direct download scripts** (provided before course)

## 9. Connectivity Test

### SSH Connection Test
```bash
# Test SSH connectivity (replace with provided details)
ssh -T username@hostname

# If successful, you should see a welcome message
```

### Internet Speed Test
Ensure you have adequate bandwidth:
- **Minimum**: 10 Mbps download
- **Recommended**: 50+ Mbps download
- **Upload**: 5+ Mbps for video calls

## 10. Backup and Recovery

### Create Backups
```bash
# Backup SSH keys
cp ~/.ssh/id_ed25519* ~/backup_location/

# Backup Git configuration
git config --list > ~/git_config_backup.txt
```

### Recovery Commands
Keep these handy in case of issues:
```bash
# Reset Git configuration
git config --global --unset-all user.name
git config --global --unset-all user.email

# Regenerate SSH keys
rm ~/.ssh/id_ed25519*
ssh-keygen -t ed25519 -C "your.email@example.com"
```

## Troubleshooting

### Common Issues

#### Git Issues
```bash
# Fix permission issues (Linux/macOS)
sudo chown -R $USER ~/.ssh/

# Reset SSH agent
eval "$(ssh-agent -k)"
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

#### Network Issues
```bash
# Test connectivity
ping google.com
curl -I https://github.com

# Check proxy settings (if behind firewall)
echo $http_proxy
echo $https_proxy
```

#### Permission Issues (Linux)
```bash
# Fix common permission problems
sudo chown -R $USER:$USER ~/
chmod 755 ~/.ssh
chmod 600 ~/.ssh/id_ed25519
chmod 644 ~/.ssh/id_ed25519.pub
```

## Getting Help

### Before the Course
- **Email**: instructors@microbial-genomics-training.org
- **Slack**: Join our pre-course channel
- **Office Hours**: Weekly setup sessions (schedule TBD)

### Documentation
- [Command Line Tutorial](https://tutorial.djangogirls.org/en/intro_to_command_line/)
- [Git Tutorial](https://git-scm.com/docs/gittutorial)
- [SSH Tutorial](https://www.ssh.com/academy/ssh)

## Final Checklist

Before the course starts, verify:

- [ ] Git installed and configured
- [ ] SSH key generated and working
- [ ] Text editor installed and functional
- [ ] Course repository cloned
- [ ] HPC access tested (if provided)
- [ ] Internet connection stable
- [ ] Backup of important configurations created

## Day 1 Preview

On the first day, we'll:
1. **Verify everyone's setup** works correctly
2. **Provide additional software** through containers
3. **Distribute datasets** and access credentials  
4. **Test HPC connections** together
5. **Troubleshoot any remaining issues**

Come prepared, but don't worry if you encounter problems – we'll sort everything out together!