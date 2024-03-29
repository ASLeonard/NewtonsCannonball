# NewtonsCannonball


## Follow up papers

 - https://www.tandfonline.com/doi/full/10.1080/24750263.2018.1558294
 - https://www.sciencedirect.com/science/article/pii/S088875431630129X
 - https://onlinelibrary.wiley.com/doi/full/10.1046/j.1365-2052.1999.00475.x?sid=nlm%3Apubmed
 - https://watermark.silverchair.com/genetics1071.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAA3AwggNsBgkqhkiG9w0BBwagggNdMIIDWQIBADCCA1IGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMYnorTZdcikh0T80GAgEQgIIDI34aMgJEanhfqtoFrFgnwvwueFhTBGmBOQJ6YeEbOFdnXezv2A7YYj1MMCJoCgTk2hyRPRam2BG3F4yY3kaU-Y7YNt9lk7jWVFgPz8oR8Qy2bW4tXnCYSicFr1dTBsiOxLp3fs-qSuAqomqrnkQwF9YC31m6o1bkKn-tOShnbVZ4pKujkdkyy5rOdYywJWOP6YY-k_61skrkoZMTrBgA6vh5HnrIzY6sqEUMPDF-SH1mQi2MHjk5zej4a9ME4_hbTSQbTa1c--DU3ZSctSkxGsGSu1VrY1EfMgHVc-e4yvjI3Riv69-hdTDxVWpSwJGOiaWrQrCfFJxae8RcLCr8csYYHnGfJL2PHVewEQ7YRKqvT6VGKfeUPUMwhe7lLBHlpjKSp9nr6Nbo-VQ42K6tKc4Hz3Zogp0cPq8h23XYmA3XGNovmmHmzVXyU0zZYxvIv1rZzrJlG2iIENgagatR-byAoqPnq3kyrluhq5TEjMjyZAT5YNF_hKmUG5tEqrQJ_h6NCQ27LEHkH8NWN-i4qJdzdJivG2Q4zChv6HYxjXcXsn1fTZiEO-Jh3780acu27Sru0W-ZhubC6mSqBaqzlnyBGk7bj91BNAb5pxMzzn8q0ow_89BFsL_3ApVFGLW93FW25djyIIHNDE9WGIhVd4O4PcY9JG_GjcNCMEa0nmf5rQThXzi6aECNQexDl7xmBuxobCcUKhwMLa9nRDpiBIucwPm3Yv3hMkh619iTrIqCilpNlkIzfUmktUWhIQgF6r_YMMAbWUiHgjGydMponUQpUWKcpNq9bLUMjfTocZ7IFehuDM-la1tfRe6dS3eeSHnpYrmZ_-MRGWeg-Q43FIQDvC4M4XyVxY2o9W6KX4DsHPA_kGmz3R0k0CqiJrysKgMVsx1O4DebP7hpkwenEiATB5_A0vM2yoScmzJDjWOX5pqPoXbk4ni7YwYLGZAHHvvWyCu6bUXYgQd8Vv09OpnkWbSeja_2wY6SYY9ZrFBHSYljpLx5OAukeOo_ATk7m_u4nDLfVisJlMdb0fle7tLWr7EM-Ce0coZ3CH7tKPs7S12m
 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC320524/pdf/nar00373-0126.pdf
 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6475130/

## Notes for the pause

Need to still filter repeats in some sensible way, and deal with circularisation or lexicographical minimization.
E.g. the two 686 bp repeats in SxG map 340.65 Mb of reads, whereas filtering out the less common copy results in the single 686 bp repeat mapping 362 Mb (+6%).
Presumably this would continue for other near-duplicate repeats or monomer-HORs.

minimap2 may also be suboptimal for the exact alignment, but first step is the filtering and decomposition.

